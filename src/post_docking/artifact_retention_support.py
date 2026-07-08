from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

from config.output_paths import output_root

_PH_RE = re.compile(r"^ph\d", re.IGNORECASE)
_VARIANTS = {"APO", "HOLO"}
_ARCHIVE_SUFFIXES = (
    ".sdf",
    ".sdf.gz",
    ".sd",
    ".pdbqt",
    ".pdbqt.gz",
    ".mol2",
    ".mol2.gz",
    ".dok",
    ".dlg",
    ".mae",
    ".maegz",
    ".ligprep_source.json",
    ".mmgbsa_pose.json",
)
_MINIMAL_EXTRA_SUFFIXES = {".tmp", ".part", ".bak", ".lock"}
_ARCHIVE_ROOT_NAME = "_artifact_archives"
_MANIFEST_FILENAME = "artifacts.manifest.json"
_ARCHIVE_FILENAME = "artifacts.tar.zst"
_INDEX_FILENAME = "archive_index.json"


@dataclass(frozen=True)
class GroupKey:
    run_id: str
    source_root: str
    pdb_id: str
    variant: str
    ph: str

    def as_dict(self) -> Dict[str, str]:
        return {
            "run_id": self.run_id,
            "source_root": self.source_root,
            "pdb_id": self.pdb_id,
            "variant": self.variant,
            "ph": self.ph,
        }

    def group_id(self) -> str:
        return (
            f"run={self.run_id}|source={self.source_root}|pdb={self.pdb_id}|"
            f"variant={self.variant}|ph={self.ph}"
        )


@dataclass(frozen=True)
class CandidateFile:
    key: GroupKey
    original_path: Path
    member_name: str
    file_type: str
    stage_dir: str
    mode: str


@dataclass
class RetentionStats:
    mode: str
    dry_run: bool
    groups_discovered: int = 0
    groups_planned: int = 0
    groups_archived: int = 0
    groups_skipped: int = 0
    groups_verified: int = 0
    groups_failed: int = 0
    files_discovered: int = 0
    files_archived: int = 0
    files_deleted: int = 0
    files_skipped: int = 0
    files_failed: int = 0
    extra_files_deleted: int = 0
    files_restored: int = 0

    def as_dict(self) -> Dict[str, Any]:
        return {
            "mode": self.mode,
            "dry_run": self.dry_run,
            "groups_discovered": self.groups_discovered,
            "groups_planned": self.groups_planned,
            "groups_archived": self.groups_archived,
            "groups_skipped": self.groups_skipped,
            "groups_verified": self.groups_verified,
            "groups_failed": self.groups_failed,
            "files_discovered": self.files_discovered,
            "files_archived": self.files_archived,
            "files_deleted": self.files_deleted,
            "files_skipped": self.files_skipped,
            "files_failed": self.files_failed,
            "extra_files_deleted": self.extra_files_deleted,
            "files_restored": self.files_restored,
        }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_inputs_safe(*, logger: Optional[logging.Logger] = None) -> Dict[str, Any]:
    try:
        from config.runtime_config import load_inputs as _load_inputs
    except Exception as exc:
        if logger is not None:
            logger.warning(
                "[artifact-retention.config] status=fallback reason=load_inputs_import_failed err=%s",
                str(exc),
            )
        return {}

    try:
        return dict(_load_inputs())
    except Exception as exc:
        if logger is not None:
            logger.warning(
                "[artifact-retention.config] status=fallback reason=load_inputs_failed err=%s",
                str(exc),
            )
        return {}


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return token or "none"


def normalize_retention_mode(raw: Any) -> str:
    token = str(raw or "").strip().lower().replace("-", "_")
    if token in {"", "rerunsafe", "rerun_safe"}:
        return "rerun_safe"
    if token in {"minimaldisk", "minimal_disk"}:
        return "minimal_disk"
    if token in {"off", "none", "false", "0", "disabled", "disable"}:
        return "off"
    return "rerun_safe"


def normalize_pdb_id(raw: Any) -> str:
    token = str(raw or "").strip().upper()
    if not token:
        return ""
    return token.replace(".PDB", "")


def normalize_variant_filter(raw: Any) -> str:
    token = str(raw or "").strip().upper()
    if not token:
        return ""
    if token in {"BASE", "NONE", "LEGACY"}:
        return "LEGACY"
    return token


def normalize_ph_filter(raw: Any) -> str:
    token = str(raw or "").strip()
    if not token:
        return ""
    if token.lower() in {"base", "none", "legacy"}:
        return "none"
    return token.lower()


def normalize_run_id(raw: Any) -> str:
    return str(raw or "").strip()


def normalize_variant_token(raw: str) -> str:
    token = normalize_variant_filter(raw)
    return token if token else "LEGACY"


def normalize_ph_token(raw: str) -> str:
    token = normalize_ph_filter(raw)
    return token if token else "none"


def parse_combo_token(raw: str) -> Tuple[str, str, str]:
    text = str(raw or "").strip()
    if not text:
        raise SystemExit("empty --combo token")
    parts = [piece.strip() for piece in text.split(":")]
    if len(parts) != 3:
        raise SystemExit(f"invalid --combo token '{text}'; expected PDB:VARIANT:PH")
    pdb_id = normalize_pdb_id(parts[0])
    if not pdb_id:
        raise SystemExit(f"invalid --combo token '{text}'; missing pdb id")
    variant = normalize_variant_token(parts[1])
    ph = normalize_ph_token(parts[2])
    return pdb_id, variant, ph


def load_combo_file(path: Path) -> List[str]:
    tokens: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens.append(stripped)
    return tokens


def resolve_combo_scope(args: argparse.Namespace) -> set[Tuple[str, str, str]]:
    raw_tokens: List[str] = [str(token) for token in list(args.combo or [])]
    combo_file = str(getattr(args, "combo_file", "") or "").strip()
    if combo_file:
        combo_path = Path(combo_file).resolve()
        if not combo_path.exists():
            raise SystemExit(f"--combo-file not found: {combo_path}")
        raw_tokens.extend(load_combo_file(combo_path))
    return {parse_combo_token(token) for token in raw_tokens}


def tooling_preflight() -> Tuple[bool, str]:
    tar_bin = shutil.which("tar")
    zstd_bin = shutil.which("zstd")
    if not tar_bin:
        return False, "tar_missing"
    if not zstd_bin:
        return False, "zstd_missing"
    try:
        proc = subprocess.run([tar_bin, "--help"], check=False, capture_output=True, text=True)
        help_text = (proc.stdout or "") + (proc.stderr or "")
    except Exception:
        return False, "tar_help_failed"
    if "--zstd" not in help_text:
        return False, "tar_missing_zstd_support"
    return True, "ok"


def match_globs(path_text: str, patterns: Sequence[str]) -> bool:
    if not patterns:
        return False
    norm = path_text.replace("\\", "/")
    return any(fnmatch.fnmatch(norm, pat) for pat in patterns)


def iter_files(root: Path) -> Iterator[Path]:
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        filenames.sort()
        for name in filenames:
            p = Path(dirpath) / name
            if p.is_file():
                yield p


def file_type(path: Path) -> str:
    lower = path.name.lower()
    for suffix in (".sdf.gz", ".pdbqt.gz", ".mol2.gz"):
        if lower.endswith(suffix):
            return suffix.lstrip(".")
    ext = path.suffix.lower().lstrip(".")
    return ext or "none"


def is_archive_target(path: Path) -> bool:
    lower = path.name.lower()
    return any(lower.endswith(sfx) for sfx in _ARCHIVE_SUFFIXES)


def infer_mode(parts: Sequence[str]) -> str:
    text = "/".join(p.lower() for p in parts)
    if "decoy" in text or "dud" in text:
        return "decoy"
    if "fda" in text:
        return "fda"
    return "unknown"


def parse_group_fields(*, run_id: str, source_root: str, rel_run_path: Path) -> GroupKey:
    parts = list(rel_run_path.parts)
    pdb_id = parts[0] if parts else "UNKNOWN"
    idx = 1
    variant = "legacy"
    if idx < len(parts) and parts[idx].upper() in _VARIANTS:
        variant = parts[idx].upper()
        idx += 1
    ph = "none"
    if idx < len(parts) and _PH_RE.match(parts[idx] or ""):
        ph = parts[idx]
    return GroupKey(run_id=run_id, source_root=source_root, pdb_id=pdb_id, variant=variant, ph=ph)


def infer_stage_dir(parts: Sequence[str]) -> str:
    stage_dir = "root"
    for piece in parts[2:-1]:
        if piece.lower().startswith("stage"):
            stage_dir = piece
            break
    if stage_dir == "root" and len(parts) >= 2:
        stage_dir = parts[-2]
    return stage_dir


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def hash_files(paths: Sequence[Path], threads: int) -> Dict[Path, str]:
    def _hash_one(path: Path) -> Tuple[Path, Optional[str]]:
        try:
            return path, hash_file(path)
        except OSError:
            return path, None

    if threads <= 1 or len(paths) <= 1:
        out_single: Dict[Path, str] = {}
        for p in paths:
            _, sha = _hash_one(p)
            if sha:
                out_single[p] = sha
        return out_single

    out: Dict[Path, str] = {}
    with ThreadPoolExecutor(max_workers=threads) as pool:
        for p, sha in pool.map(_hash_one, paths):
            if sha:
                out[p] = sha
    return out


def archive_group_path(archive_root: Path, key: GroupKey) -> Path:
    group_dir = (
        archive_root
        / safe_token(key.source_root)
        / safe_token(key.pdb_id)
        / safe_token(key.variant)
        / safe_token(key.ph)
    )
    return group_dir / _ARCHIVE_FILENAME


def write_tar(*, entries: Sequence[Dict[str, Any]], path_by_member: Mapping[str, Path], tar_path: Path) -> None:
    tar_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, mode="w", format=tarfile.PAX_FORMAT) as tf:
        for entry in sorted(entries, key=lambda x: str(x["member_name"])):
            member_name = str(entry["member_name"])
            src = path_by_member[member_name]
            info = tf.gettarinfo(str(src), arcname=member_name)
            info.uid = 0
            info.gid = 0
            info.uname = ""
            info.gname = ""
            with src.open("rb") as fh:
                tf.addfile(info, fh)


def compress_tar_to_zst(tar_path: Path, archive_path: Path, threads: int) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_archive = archive_path.with_suffix(".tar.zst.tmp")
    cmd = ["zstd", "-q", "-f", f"-T{max(1, int(threads))}", "--rm", str(tar_path), "-o", str(tmp_archive)]
    subprocess.run(cmd, check=True)
    os.replace(tmp_archive, archive_path)


def verify_archive(archive_path: Path, entries: Sequence[Dict[str, Any]]) -> Tuple[bool, List[str]]:
    errors: List[str] = []
    if not archive_path.exists():
        return False, [f"archive_missing:{archive_path}"]
    with tempfile.TemporaryDirectory(prefix="artifact-retention-verify-") as tmp:
        tmp_root = Path(tmp)
        proc = subprocess.run(["tar", "--zstd", "-xf", str(archive_path), "-C", str(tmp_root)], check=False)
        if proc.returncode != 0:
            return False, [f"extract_failed:returncode={proc.returncode}"]
        for entry in entries:
            member_name = str(entry["member_name"])
            expected_size = int(entry["size_bytes"])
            expected_sha = str(entry["sha256"])
            restored = tmp_root / member_name
            if not restored.exists():
                errors.append(f"missing_member:{member_name}")
                continue
            if not restored.is_file():
                errors.append(f"not_file_member:{member_name}")
                continue
            size_now = restored.stat().st_size
            if size_now != expected_size:
                errors.append(f"size_mismatch:{member_name}:expected={expected_size}:found={size_now}")
                continue
            sha_now = hash_file(restored)
            if sha_now != expected_sha:
                errors.append(f"sha_mismatch:{member_name}")
    return len(errors) == 0, errors


def load_manifest(manifest_path: Path) -> Dict[str, Any]:
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def delete_files(paths: Sequence[Path]) -> Tuple[int, int]:
    deleted = 0
    failed = 0
    for p in paths:
        try:
            if p.exists():
                p.unlink()
                deleted += 1
        except Exception:
            failed += 1
    return deleted, failed


def discover_candidates(*, run_id: str, run_docked: Path, run_post_docked: Path, pdb_id_filter: str, variant_filter: str, ph_filter: str, include_globs: Sequence[str], exclude_globs: Sequence[str], logger: logging.Logger) -> Dict[GroupKey, List[CandidateFile]]:
    groups: Dict[GroupKey, List[CandidateFile]] = {}
    for source_root, root in (("docked", run_docked), ("post_docked", run_post_docked)):
        if not root.is_dir():
            logger.info("[artifact-retention.discover] action=skip source=%s reason=missing_root path=%s", source_root, root)
            continue
        for abs_path in iter_files(root):
            rel = abs_path.relative_to(root)
            member_name = f"{source_root}/{run_id}/{rel.as_posix()}"
            if f"/{_ARCHIVE_ROOT_NAME}/" in abs_path.as_posix() or abs_path.name == _MANIFEST_FILENAME:
                continue
            if abs_path.suffix.lower() == ".csv" or not is_archive_target(abs_path):
                continue
            if include_globs and not match_globs(member_name, include_globs):
                continue
            if match_globs(member_name, exclude_globs):
                continue
            key = parse_group_fields(run_id=run_id, source_root=source_root, rel_run_path=rel)
            if pdb_id_filter and key.pdb_id.upper() != pdb_id_filter:
                continue
            if variant_filter and key.variant.upper() != variant_filter:
                continue
            if ph_filter and str(key.ph).strip().lower() != ph_filter:
                continue
            item = CandidateFile(
                key=key,
                original_path=abs_path,
                member_name=member_name,
                file_type=file_type(abs_path),
                stage_dir=infer_stage_dir(rel.parts),
                mode=infer_mode(rel.parts),
            )
            groups.setdefault(key, []).append(item)
    for key in list(groups.keys()):
        groups[key] = sorted(groups[key], key=lambda x: x.member_name)
    return groups


def collect_extra_minimal_disk_files(*, run_id: str, run_docked: Path, run_post_docked: Path, pdb_id_filter: str, variant_filter: str, ph_filter: str, include_globs: Sequence[str], exclude_globs: Sequence[str]) -> List[Path]:
    out: List[Path] = []
    for source_root, root in (("docked", run_docked), ("post_docked", run_post_docked)):
        if not root.is_dir():
            continue
        for abs_path in iter_files(root):
            lower_name = abs_path.name.lower()
            if abs_path.suffix.lower() == ".csv" or f"/{_ARCHIVE_ROOT_NAME}/" in abs_path.as_posix():
                continue
            rel = abs_path.relative_to(root)
            key = parse_group_fields(run_id=run_id, source_root=source_root, rel_run_path=rel)
            if pdb_id_filter and key.pdb_id.upper() != pdb_id_filter:
                continue
            if variant_filter and key.variant.upper() != variant_filter:
                continue
            if ph_filter and str(key.ph).strip().lower() != ph_filter:
                continue
            if abs_path.suffix.lower() not in _MINIMAL_EXTRA_SUFFIXES:
                continue
            member_name = f"{source_root}/{run_id}/{rel.as_posix()}"
            if include_globs and not match_globs(member_name, include_globs):
                continue
            if match_globs(member_name, exclude_globs) or lower_name.endswith('.csv'):
                continue
            out.append(abs_path)
    out.sort()
    return out


def load_run_index(index_path: Path) -> Dict[str, Any]:
    if not index_path.exists():
        return {"schema_version": 2, "groups": []}
    try:
        data = json.loads(index_path.read_text(encoding="utf-8"))
    except Exception:
        return {"schema_version": 2, "groups": []}
    if not isinstance(data, dict):
        return {"schema_version": 2, "groups": []}
    if not isinstance(data.get("groups"), list):
        data["groups"] = []
    return data


def write_run_index(index_path: Path, *, run_id: str, mode: str, groups_map: Mapping[str, Dict[str, Any]]) -> None:
    groups = [groups_map[k] for k in sorted(groups_map.keys())]
    payload = {
        "schema_version": 2,
        "generated_at": utc_now(),
        "run_id": run_id,
        "mode": mode,
        "group_count": len(groups),
        "groups": groups,
    }
    index_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = index_path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=False), encoding='utf-8')
    os.replace(tmp, index_path)


def build_entries(files: Sequence[CandidateFile], archive_path: Path, threads: int) -> List[Dict[str, Any]]:
    hashes = hash_files([f.original_path for f in files], threads=threads)
    entries: List[Dict[str, Any]] = []
    for f in files:
        sha = hashes.get(f.original_path)
        if not sha:
            continue
        try:
            st = f.original_path.stat()
        except OSError:
            continue
        entries.append({
            "run_id": f.key.run_id,
            "pdb_id": f.key.pdb_id,
            "variant": f.key.variant,
            "ph": f.key.ph,
            "stage_dir": f.stage_dir,
            "mode": f.mode,
            "original_path": str(f.original_path),
            "archive_path": str(archive_path),
            "member_name": f.member_name,
            "size_bytes": int(st.st_size),
            "mtime": float(st.st_mtime),
            "sha256": sha,
            "file_type": f.file_type,
        })
    entries.sort(key=lambda x: str(x["member_name"]))
    return entries


def resolve_run_roots(args: argparse.Namespace, *, repo_root: Path, cfg: Mapping[str, Any], run_id: str) -> Tuple[Path, Path, Path, Path]:
    overall = Path(str(cfg.get("OVERALL_DIR", repo_root))).resolve()
    docked_base = (
        Path(args.docked_root).resolve()
        if args.docked_root
        else Path(str(cfg.get("DOCKED_DIR", output_root(overall, "docked")))).resolve()
    )
    post_base = (
        Path(args.post_docked_root).resolve()
        if args.post_docked_root
        else Path(
            str(cfg.get("POST_DOCKED_DIR", output_root(overall, "post_docked")))
        ).resolve()
    )
    run_docked = docked_base if docked_base.name == run_id else docked_base / run_id
    run_post = post_base if post_base.name == run_id else post_base / run_id
    return docked_base, post_base, run_docked, run_post


def resolve_archive_path(raw_path: Any, base_dir: Path) -> Path:
    p = Path(str(raw_path or "").strip()).expanduser()
    return ((base_dir / p).resolve() if not p.is_absolute() else p.resolve())


def build_group_record(*, key: GroupKey, archive_path: Path, entries: Sequence[Dict[str, Any]], verified: bool, verify_errors: Sequence[str]) -> Dict[str, Any]:
    return {
        "group_id": key.group_id(),
        "group": key.as_dict(),
        "archive_path": str(archive_path),
        "entry_count": len(entries),
        "verified": bool(verified),
        "verify_errors": list(verify_errors),
        "entries": list(entries),
    }


def is_restore_target_safe(target: Path, *, allowed_root_a: Path, allowed_root_b: Path) -> bool:
    try:
        resolved = target.resolve()
    except Exception:
        return False
    try:
        resolved.relative_to(allowed_root_a.resolve())
        return True
    except Exception:
        pass
    try:
        resolved.relative_to(allowed_root_b.resolve())
        return True
    except Exception:
        pass
    return False


def restore_archive(archive_path: Path, *, payload_label: str, payload: Mapping[str, Any], allowed_root_a: Path, allowed_root_b: Path, logger: logging.Logger) -> Tuple[int, int]:
    entries = payload.get("entries")
    if not isinstance(entries, list) or not entries:
        raise SystemExit(f"restore payload has no entries: {payload_label}")
    restored = 0
    failed = 0
    with tempfile.TemporaryDirectory(prefix="artifact-retention-restore-") as tmp:
        tmp_root = Path(tmp)
        proc = subprocess.run(["tar", "--zstd", "-xf", str(archive_path), "-C", str(tmp_root)], check=False)
        if proc.returncode != 0:
            raise SystemExit(f"failed to extract archive: {archive_path} (rc={proc.returncode})")
        for entry in entries:
            if not isinstance(entry, dict):
                failed += 1
                continue
            member_name = str(entry.get("member_name", "")).strip()
            original_path = str(entry.get("original_path", "")).strip()
            if not member_name or not original_path:
                failed += 1
                continue
            src = tmp_root / member_name
            dst = Path(original_path).resolve()
            if not src.exists() or not src.is_file():
                logger.error("[artifact-retention.restore] member_missing member=%s archive=%s", member_name, archive_path)
                failed += 1
                continue
            if not is_restore_target_safe(dst, allowed_root_a=allowed_root_a, allowed_root_b=allowed_root_b):
                logger.error("[artifact-retention.restore] target_unsafe path=%s archive=%s", dst, archive_path)
                failed += 1
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            try:
                shutil.copy2(src, dst)
                restored += 1
            except Exception:
                logger.error("[artifact-retention.restore] copy_failed src=%s dst=%s", src, dst, exc_info=True)
                failed += 1
    return restored, failed


def resolve_archive_root_for_run(args: argparse.Namespace, run_id: str) -> Path:
    if args.docked_root:
        docked_root = Path(args.docked_root).resolve()
        run_docked = docked_root if docked_root.name == run_id else docked_root / run_id
        return run_docked / _ARCHIVE_ROOT_NAME
    repo_root = Path(args.repo_root).resolve() if args.repo_root else Path(__file__).resolve().parents[2]
    return (output_root(repo_root, "docked") / run_id / _ARCHIVE_ROOT_NAME).resolve()


def infer_run_id_from_archive_path(archive_path: Path) -> Optional[str]:
    parts = list(archive_path.resolve().parts)
    for idx, piece in enumerate(parts):
        if piece == _ARCHIVE_ROOT_NAME and idx >= 2:
            return parts[idx - 1]
    return None


def resolve_restore_payload(args: argparse.Namespace, *, run_id: Optional[str]) -> Tuple[Path, Dict[str, Any], str]:
    if args.restore_group:
        if not run_id:
            raise SystemExit("--restore-group requires --run-id")
        archive_root = resolve_archive_root_for_run(args, run_id)
        index_path = archive_root / _INDEX_FILENAME
        if not index_path.exists():
            raise SystemExit(f"archive index not found: {index_path}")
        index_payload = load_run_index(index_path)
        for item in index_payload.get("groups", []):
            if str(item.get("group_id", "")).strip() != str(args.restore_group).strip():
                continue
            archive_path = resolve_archive_path(item.get("archive_path", ""), archive_root)
            entries = item.get("entries")
            if isinstance(entries, list) and entries:
                payload = {
                    "group": item.get("group") or {},
                    "group_id": item.get("group_id"),
                    "archive_path": str(archive_path),
                    "entries": entries,
                }
                return archive_path, payload, f"{index_path}#{args.restore_group}"
            manifest_text = str(item.get("manifest_path", "")).strip()
            if manifest_text:
                manifest_path = Path(manifest_text).resolve()
                if not manifest_path.exists():
                    raise SystemExit(f"legacy manifest not found: {manifest_path}")
                return archive_path, load_manifest(manifest_path), str(manifest_path)
            raise SystemExit(f"group has no restore metadata: {args.restore_group}")
        raise SystemExit(f"group not found in index: {args.restore_group}")
    if not args.archive:
        raise SystemExit("--restore requires --archive <path> or --restore-group <group-id> with --run-id")
    archive_path = Path(args.archive).resolve()
    if args.manifest:
        manifest_path = Path(args.manifest).resolve()
        if not manifest_path.exists():
            raise SystemExit(f"manifest not found: {manifest_path}")
        return archive_path, load_manifest(manifest_path), str(manifest_path)
    candidate_run_id = run_id or infer_run_id_from_archive_path(archive_path)
    if candidate_run_id:
        archive_root = resolve_archive_root_for_run(args, candidate_run_id)
        index_path = archive_root / _INDEX_FILENAME
        if index_path.exists():
            index_payload = load_run_index(index_path)
            archive_norm = archive_path.resolve()
            for item in index_payload.get("groups", []):
                item_archive = resolve_archive_path(item.get("archive_path", ""), archive_root)
                if item_archive != archive_norm:
                    continue
                entries = item.get("entries")
                if isinstance(entries, list) and entries:
                    payload = {
                        "group": item.get("group") or {},
                        "group_id": item.get("group_id"),
                        "archive_path": str(archive_path),
                        "entries": entries,
                    }
                    return archive_path, payload, f"{index_path}#{item.get('group_id')}"
                manifest_text = str(item.get("manifest_path", "")).strip()
                if manifest_text:
                    manifest_path = Path(manifest_text).resolve()
                    if manifest_path.exists():
                        return archive_path, load_manifest(manifest_path), str(manifest_path)
    legacy_manifest = archive_path.with_name(_MANIFEST_FILENAME)
    if legacy_manifest.exists():
        return archive_path, load_manifest(legacy_manifest), str(legacy_manifest)
    raise SystemExit(
        f"unable to resolve restore metadata for archive: {archive_path}; use --manifest or --restore-group with --run-id"
    )


def filter_restore_entries(entries: Sequence[Dict[str, Any]], *, stage: Optional[str], mode: Optional[str], path_globs: Sequence[str]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    wanted_stage = str(stage or "").strip()
    wanted_mode = str(mode or "").strip().lower()
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        if wanted_stage and str(entry.get("stage_dir", "")) != wanted_stage:
            continue
        entry_mode = str(entry.get("mode", "")).strip().lower()
        if wanted_mode and entry_mode != wanted_mode:
            continue
        if path_globs:
            member = str(entry.get("member_name", ""))
            original = str(entry.get("original_path", ""))
            if not any(fnmatch.fnmatch(member, pat) or fnmatch.fnmatch(original, pat) for pat in path_globs):
                continue
        out.append(entry)
    return out
