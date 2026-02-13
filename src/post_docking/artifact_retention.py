# -*- coding: utf-8 -*-
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

from input_and_export_functions import load_inputs


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


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")


def _to_bool(raw: Any, default: bool = False) -> bool:
    if raw is None:
        return default
    if isinstance(raw, bool):
        return raw
    token = str(raw).strip().lower()
    if token in {"1", "true", "yes", "on", "y"}:
        return True
    if token in {"0", "false", "no", "off", "n"}:
        return False
    return default


def _safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return token or "none"


def _normalize_retention_mode(raw: Any) -> str:
    token = str(raw or "").strip().lower().replace("-", "_")
    if token in {"", "rerunsafe", "rerun_safe"}:
        return "rerun_safe"
    if token in {"minimaldisk", "minimal_disk"}:
        return "minimal_disk"
    if token in {"off", "none", "false", "0", "disabled", "disable"}:
        return "off"
    return "rerun_safe"


def _normalize_pdb_id(raw: Any) -> str:
    token = str(raw or "").strip().upper()
    if not token:
        return ""
    return token.replace(".PDB", "")


def _normalize_run_id(raw: Any) -> str:
    return str(raw or "").strip()


def _tooling_preflight() -> Tuple[bool, str]:
    tar_bin = shutil.which("tar")
    zstd_bin = shutil.which("zstd")
    if not tar_bin:
        return False, "tar_missing"
    if not zstd_bin:
        return False, "zstd_missing"
    try:
        proc = subprocess.run(
            [tar_bin, "--help"],
            check=False,
            capture_output=True,
            text=True,
        )
        help_text = (proc.stdout or "") + (proc.stderr or "")
    except Exception:
        return False, "tar_help_failed"
    if "--zstd" not in help_text:
        return False, "tar_missing_zstd_support"
    return True, "ok"


def _match_globs(path_text: str, patterns: Sequence[str]) -> bool:
    if not patterns:
        return False
    norm = path_text.replace("\\", "/")
    for pat in patterns:
        if fnmatch.fnmatch(norm, pat):
            return True
    return False


def _iter_files(root: Path) -> Iterator[Path]:
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        filenames.sort()
        for name in filenames:
            p = Path(dirpath) / name
            if p.is_file():
                yield p


def _file_type(path: Path) -> str:
    lower = path.name.lower()
    for suffix in (".sdf.gz", ".pdbqt.gz", ".mol2.gz"):
        if lower.endswith(suffix):
            return suffix.lstrip(".")
    ext = path.suffix.lower().lstrip(".")
    return ext or "none"


def _is_archive_target(path: Path) -> bool:
    lower = path.name.lower()
    return any(lower.endswith(sfx) for sfx in _ARCHIVE_SUFFIXES)


def _infer_mode(parts: Sequence[str]) -> str:
    text = "/".join(p.lower() for p in parts)
    if "decoy" in text or "dud" in text:
        return "decoy"
    if "fda" in text:
        return "fda"
    return "unknown"


def _parse_group_fields(
    *,
    run_id: str,
    source_root: str,
    rel_run_path: Path,
) -> GroupKey:
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
        idx += 1

    return GroupKey(
        run_id=run_id,
        source_root=source_root,
        pdb_id=pdb_id,
        variant=variant,
        ph=ph,
    )


def _infer_stage_dir(parts: Sequence[str]) -> str:
    stage_dir = "root"
    for piece in parts[2:-1]:
        if piece.lower().startswith("stage"):
            stage_dir = piece
            break
    if stage_dir == "root" and len(parts) >= 2:
        stage_dir = parts[-2]
    return stage_dir


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _hash_files(paths: Sequence[Path], threads: int) -> Dict[Path, str]:
    def _hash_one(path: Path) -> Tuple[Path, Optional[str]]:
        try:
            return path, _hash_file(path)
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


def _archive_group_path(archive_root: Path, key: GroupKey) -> Path:
    group_dir = (
        archive_root
        / _safe_token(key.source_root)
        / _safe_token(key.pdb_id)
        / _safe_token(key.variant)
        / _safe_token(key.ph)
    )
    return group_dir / _ARCHIVE_FILENAME


def _write_tar(
    *,
    entries: Sequence[Dict[str, Any]],
    path_by_member: Mapping[str, Path],
    tar_path: Path,
) -> None:
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


def _compress_tar_to_zst(tar_path: Path, archive_path: Path, threads: int) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_archive = archive_path.with_suffix(".tar.zst.tmp")
    cmd = [
        "zstd",
        "-q",
        "-f",
        f"-T{max(1, int(threads))}",
        "--rm",
        str(tar_path),
        "-o",
        str(tmp_archive),
    ]
    subprocess.run(cmd, check=True)
    os.replace(tmp_archive, archive_path)


def _verify_archive(
    archive_path: Path,
    entries: Sequence[Dict[str, Any]],
) -> Tuple[bool, List[str]]:
    errors: List[str] = []
    if not archive_path.exists():
        return False, [f"archive_missing:{archive_path}"]

    with tempfile.TemporaryDirectory(prefix="artifact-retention-verify-") as tmp:
        tmp_root = Path(tmp)
        cmd = ["tar", "--zstd", "-xf", str(archive_path), "-C", str(tmp_root)]
        proc = subprocess.run(cmd, check=False)
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
                errors.append(
                    f"size_mismatch:{member_name}:expected={expected_size}:found={size_now}"
                )
                continue
            sha_now = _hash_file(restored)
            if sha_now != expected_sha:
                errors.append(f"sha_mismatch:{member_name}")

    return (len(errors) == 0), errors


def _load_manifest(manifest_path: Path) -> Dict[str, Any]:
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def _delete_files(paths: Sequence[Path]) -> Tuple[int, int]:
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


def _discover_candidates(
    *,
    run_id: str,
    run_docked: Path,
    run_post_docked: Path,
    pdb_id_filter: str,
    include_globs: Sequence[str],
    exclude_globs: Sequence[str],
    logger: logging.Logger,
) -> Dict[GroupKey, List[CandidateFile]]:
    groups: Dict[GroupKey, List[CandidateFile]] = {}
    roots = [("docked", run_docked), ("post_docked", run_post_docked)]

    for source_root, root in roots:
        if not root.is_dir():
            logger.info(
                "[artifact-retention.discover] action=skip source=%s reason=missing_root path=%s",
                source_root,
                root,
            )
            continue

        for abs_path in _iter_files(root):
            rel = abs_path.relative_to(root)
            member_name = f"{source_root}/{run_id}/{rel.as_posix()}"

            if f"/{_ARCHIVE_ROOT_NAME}/" in abs_path.as_posix():
                continue
            if abs_path.name == _MANIFEST_FILENAME:
                continue
            if abs_path.suffix.lower() == ".csv":
                continue
            if not _is_archive_target(abs_path):
                continue
            if include_globs and not _match_globs(member_name, include_globs):
                continue
            if _match_globs(member_name, exclude_globs):
                continue

            key = _parse_group_fields(
                run_id=run_id,
                source_root=source_root,
                rel_run_path=rel,
            )
            if pdb_id_filter and key.pdb_id.upper() != pdb_id_filter:
                continue
            item = CandidateFile(
                key=key,
                original_path=abs_path,
                member_name=member_name,
                file_type=_file_type(abs_path),
                stage_dir=_infer_stage_dir(rel.parts),
                mode=_infer_mode(rel.parts),
            )
            groups.setdefault(key, []).append(item)

    for key in list(groups.keys()):
        groups[key] = sorted(groups[key], key=lambda x: x.member_name)
    return groups


def _collect_extra_minimal_disk_files(
    *,
    run_id: str,
    run_docked: Path,
    run_post_docked: Path,
    pdb_id_filter: str,
    include_globs: Sequence[str],
    exclude_globs: Sequence[str],
) -> List[Path]:
    out: List[Path] = []
    for source_root, root in (("docked", run_docked), ("post_docked", run_post_docked)):
        if not root.is_dir():
            continue
        for abs_path in _iter_files(root):
            lower_name = abs_path.name.lower()
            if abs_path.suffix.lower() == ".csv":
                continue
            if f"/{_ARCHIVE_ROOT_NAME}/" in abs_path.as_posix():
                continue
            rel = abs_path.relative_to(root)
            rel_parts = rel.parts
            rel_pdb = rel_parts[0].upper() if rel_parts else ""
            if pdb_id_filter and rel_pdb != pdb_id_filter:
                continue
            if abs_path.suffix.lower() not in _MINIMAL_EXTRA_SUFFIXES:
                continue
            member_name = f"{source_root}/{run_id}/{rel.as_posix()}"
            if include_globs and not _match_globs(member_name, include_globs):
                continue
            if _match_globs(member_name, exclude_globs):
                continue
            if lower_name.endswith(".csv"):
                continue
            out.append(abs_path)
    out.sort()
    return out


def _load_run_index(index_path: Path) -> Dict[str, Any]:
    if not index_path.exists():
        return {"schema_version": 2, "groups": []}
    try:
        data = json.loads(index_path.read_text(encoding="utf-8"))
    except Exception:
        return {"schema_version": 2, "groups": []}
    if not isinstance(data, dict):
        return {"schema_version": 2, "groups": []}
    groups = data.get("groups")
    if not isinstance(groups, list):
        data["groups"] = []
    return data


def _write_run_index(
    index_path: Path,
    *,
    run_id: str,
    mode: str,
    groups_map: Mapping[str, Dict[str, Any]],
) -> None:
    groups = [groups_map[k] for k in sorted(groups_map.keys())]
    payload = {
        "schema_version": 2,
        "generated_at": _utc_now(),
        "run_id": run_id,
        "mode": mode,
        "group_count": len(groups),
        "groups": groups,
    }
    index_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = index_path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=False), encoding="utf-8")
    os.replace(tmp, index_path)


def _build_entries(
    files: Sequence[CandidateFile],
    archive_path: Path,
    threads: int,
) -> List[Dict[str, Any]]:
    paths = [f.original_path for f in files]
    hashes = _hash_files(paths, threads=threads)
    entries: List[Dict[str, Any]] = []
    for f in files:
        sha = hashes.get(f.original_path)
        if not sha:
            continue
        try:
            st = f.original_path.stat()
        except OSError:
            continue
        entry = {
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
        }
        entries.append(entry)
    entries.sort(key=lambda x: str(x["member_name"]))
    return entries


def _resolve_run_roots(
    args: argparse.Namespace,
    *,
    repo_root: Path,
    cfg: Mapping[str, Any],
    run_id: str,
) -> Tuple[Path, Path, Path, Path]:
    overall = Path(str(cfg.get("OVERALL_DIR", repo_root))).resolve()
    docked_base = (
        Path(args.docked_root).resolve()
        if args.docked_root
        else Path(str(cfg.get("DOCKED_DIR", overall / "docked"))).resolve()
    )
    post_base = (
        Path(args.post_docked_root).resolve()
        if args.post_docked_root
        else Path(str(cfg.get("POST_DOCKED_DIR", overall / "post_docked"))).resolve()
    )
    run_docked = docked_base if docked_base.name == run_id else docked_base / run_id
    run_post = post_base if post_base.name == run_id else post_base / run_id
    return docked_base, post_base, run_docked, run_post


def _resolve_archive_path(raw_path: Any, base_dir: Path) -> Path:
    p = Path(str(raw_path or "").strip()).expanduser()
    if not p.is_absolute():
        p = (base_dir / p).resolve()
    else:
        p = p.resolve()
    return p


def _build_group_record(
    *,
    key: GroupKey,
    archive_path: Path,
    entries: Sequence[Dict[str, Any]],
    verified: bool,
    verify_errors: Sequence[str],
) -> Dict[str, Any]:
    return {
        "group_id": key.group_id(),
        "group": key.as_dict(),
        "archive_path": str(archive_path),
        "entry_count": len(entries),
        "verified": bool(verified),
        "verify_errors": list(verify_errors),
        "entries": list(entries),
    }


def _run_retention(args: argparse.Namespace) -> int:
    logger = logging.getLogger("artifact-retention")
    run_id = str(args.run_id or "").strip()
    if not run_id:
        raise SystemExit("--run-id is required unless --restore with --archive is used")
    mode = _normalize_retention_mode(args.mode)
    if mode == "off":
        logger.info("[artifact-retention.skip] run_id=%s reason=mode_off", run_id)
        return 0
    pdb_id_filter = _normalize_pdb_id(args.pdb_id)

    repo_root = (
        Path(args.repo_root).resolve()
        if args.repo_root
        else Path(__file__).resolve().parents[2]
    )
    cfg: Dict[str, Any] = {}
    try:
        cfg = dict(load_inputs())
    except Exception:
        cfg = {}
    _docked_base, _post_base, run_docked, run_post = _resolve_run_roots(
        args,
        repo_root=repo_root,
        cfg=cfg,
        run_id=run_id,
    )
    archive_root = run_docked / _ARCHIVE_ROOT_NAME
    archive_root.mkdir(parents=True, exist_ok=True)

    include_globs = list(args.include_glob or [])
    exclude_globs = list(args.exclude_glob or [])
    exclude_globs.extend(
        [
            f"docked/{run_id}/{_ARCHIVE_ROOT_NAME}/*",
            f"post_docked/{run_id}/{_ARCHIVE_ROOT_NAME}/*",
        ]
    )

    stats = RetentionStats(mode=mode, dry_run=bool(args.dry_run))
    logger.info(
        "[artifact-retention.start] run_id=%s pdb_id=%s mode=%s dry_run=%s verify_only=%s overwrite=%s docked_root=%s post_docked_root=%s",
        run_id,
        pdb_id_filter or "*",
        mode,
        str(bool(args.dry_run)).lower(),
        str(bool(args.verify_only)).lower(),
        str(bool(args.overwrite)).lower(),
        run_docked,
        run_post,
    )

    preflight_ok, preflight_reason = _tooling_preflight()
    if not preflight_ok:
        logger.warning(
            "[artifact-retention.skip] run_id=%s reason=tooling_preflight_failed detail=%s",
            run_id,
            preflight_reason,
        )
        logger.info("[artifact-retention.summary] %s", json.dumps(stats.as_dict(), sort_keys=True))
        return 0

    index_path = archive_root / _INDEX_FILENAME
    index_payload = _load_run_index(index_path)
    index_groups_raw = index_payload.get("groups") or []
    index_groups = index_groups_raw if isinstance(index_groups_raw, list) else []
    groups_map: Dict[str, Dict[str, Any]] = {}
    for item in index_groups:
        if not isinstance(item, dict):
            continue
        group_id = str(item.get("group_id", "")).strip()
        if not group_id:
            continue
        groups_map[group_id] = dict(item)

    if args.verify_only:
        if index_groups:
            for item in index_groups:
                if not isinstance(item, dict):
                    continue
                group_id = str(item.get("group_id", "")).strip() or "(missing_group_id)"
                archive_path = _resolve_archive_path(item.get("archive_path", ""), archive_root)
                entries = item.get("entries")
                if not isinstance(entries, list):
                    entries = []
                if not entries:
                    manifest_text = str(item.get("manifest_path", "")).strip()
                    if manifest_text:
                        try:
                            manifest_payload = _load_manifest(Path(manifest_text).resolve())
                            entries = manifest_payload.get("entries") or []
                        except Exception:
                            stats.groups_failed += 1
                            logger.error(
                                "[artifact-retention.verify] status=failed group_id=%s reason=legacy_manifest_unreadable",
                                group_id,
                            )
                            continue
                if not entries:
                    stats.groups_skipped += 1
                    logger.warning(
                        "[artifact-retention.verify] status=skip group_id=%s reason=missing_entries",
                        group_id,
                    )
                    continue
                ok, verify_errors = _verify_archive(archive_path, entries)
                if ok:
                    stats.groups_verified += 1
                    logger.info(
                        "[artifact-retention.verify] status=ok group_id=%s archive=%s entries=%d",
                        group_id,
                        archive_path,
                        len(entries),
                    )
                else:
                    stats.groups_failed += 1
                    logger.error(
                        "[artifact-retention.verify] status=failed group_id=%s archive=%s errors=%s",
                        group_id,
                        archive_path,
                        ";".join(verify_errors[:10]),
                    )
        else:
            manifests = sorted(archive_root.rglob(_MANIFEST_FILENAME))
            if not manifests:
                logger.warning(
                    "[artifact-retention.verify] run_id=%s status=empty archive_root=%s",
                    run_id,
                    archive_root,
                )
            for manifest_path in manifests:
                try:
                    payload = _load_manifest(manifest_path)
                except Exception:
                    stats.groups_failed += 1
                    continue
                archive_path = Path(str(payload.get("archive_path", "")))
                entries = payload.get("entries") or []
                if not archive_path.is_absolute():
                    archive_path = (manifest_path.parent / archive_path).resolve()
                ok, verify_errors = _verify_archive(archive_path, entries)
                if ok:
                    stats.groups_verified += 1
                    logger.info(
                        "[artifact-retention.verify] status=ok manifest=%s archive=%s entries=%d",
                        manifest_path,
                        archive_path,
                        len(entries),
                    )
                else:
                    stats.groups_failed += 1
                    logger.error(
                        "[artifact-retention.verify] status=failed manifest=%s archive=%s errors=%s",
                        manifest_path,
                        archive_path,
                        ";".join(verify_errors[:10]),
                    )
        logger.info("[artifact-retention.summary] %s", json.dumps(stats.as_dict(), sort_keys=True))
        return 0 if stats.groups_failed == 0 else 1

    groups = _discover_candidates(
        run_id=run_id,
        run_docked=run_docked,
        run_post_docked=run_post,
        pdb_id_filter=pdb_id_filter,
        include_globs=include_globs,
        exclude_globs=exclude_globs,
        logger=logger,
    )
    stats.groups_discovered = len(groups)
    stats.files_discovered = sum(len(v) for v in groups.values())

    for key in sorted(groups.keys(), key=lambda k: k.group_id()):
        files = groups[key]
        archive_path = _archive_group_path(archive_root, key)
        stats.groups_planned += 1
        discovered_members = {item.member_name for item in files}
        existing = groups_map.get(key.group_id()) or {}
        existing_entries_raw = existing.get("entries")
        existing_entries = (
            existing_entries_raw if isinstance(existing_entries_raw, list) else []
        )
        existing_members = {
            str(entry.get("member_name", ""))
            for entry in existing_entries
            if isinstance(entry, dict)
        }
        existing_archive = _resolve_archive_path(
            existing.get("archive_path", archive_path),
            archive_root,
        )

        if not args.overwrite and bool(existing.get("verified", False)) and existing_entries:
            manifest_matches = (
                len(existing_entries) == len(files)
                and existing_members == discovered_members
            )
            if not existing_archive.exists():
                logger.warning(
                    "[artifact-retention.group] action=rebuild group_id=%s reason=archive_missing",
                    key.group_id(),
                )
            elif manifest_matches:
                stats.groups_skipped += 1
                stats.files_skipped += len(files)
                logger.info(
                    "[artifact-retention.group] action=skip group_id=%s reason=already_archived_verified files=%d",
                    key.group_id(),
                    len(files),
                )
                continue
            else:
                logger.warning(
                    "[artifact-retention.group] action=rebuild group_id=%s reason=manifest_mismatch discovered=%d manifest=%d",
                    key.group_id(),
                    len(files),
                    len(existing_entries),
                )

        if args.dry_run:
            logger.info(
                "[artifact-retention.group] action=dry_run_plan group_id=%s files=%d archive=%s index=%s",
                key.group_id(),
                len(files),
                archive_path,
                index_path,
            )
            continue

        entries = _build_entries(
            files,
            archive_path=archive_path,
            threads=max(1, int(args.threads)),
        )
        if len(entries) != len(files):
            vanished = len(files) - len(entries)
            if vanished > 0:
                stats.files_skipped += vanished
                logger.warning(
                    "[artifact-retention.group] action=partial group_id=%s reason=files_missing_during_hash_or_stat discovered=%d stable=%d",
                    key.group_id(),
                    len(files),
                    len(entries),
                )
        if not entries:
            stats.groups_skipped += 1
            logger.info(
                "[artifact-retention.group] action=skip group_id=%s reason=no_stable_files",
                key.group_id(),
            )
            continue
        stable_members = {str(entry["member_name"]) for entry in entries}
        path_by_member = {
            item.member_name: item.original_path
            for item in files
            if item.member_name in stable_members
        }

        archive_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_tar = archive_path.with_suffix(".tar.tmp")
        if tmp_tar.exists():
            tmp_tar.unlink()
        if archive_path.exists() and args.overwrite:
            archive_path.unlink()

        try:
            _write_tar(entries=entries, path_by_member=path_by_member, tar_path=tmp_tar)
            _compress_tar_to_zst(tmp_tar, archive_path=archive_path, threads=max(1, int(args.threads)))
            ok, verify_errors = _verify_archive(archive_path, entries)
            groups_map[key.group_id()] = _build_group_record(
                key=key,
                archive_path=archive_path,
                entries=entries,
                verified=ok,
                verify_errors=verify_errors,
            )
            if not ok:
                stats.groups_failed += 1
                stats.files_failed += len(entries)
                logger.error(
                    "[artifact-retention.group] status=failed group_id=%s reason=verify_failed errors=%s",
                    key.group_id(),
                    ";".join(verify_errors[:10]),
                )
                continue

            stats.groups_archived += 1
            stats.groups_verified += 1
            stats.files_archived += len(entries)
            to_delete = list(path_by_member.values())
            deleted, failed = _delete_files(to_delete)
            stats.files_deleted += deleted
            stats.files_failed += failed
            logger.info(
                "[artifact-retention.group] status=ok group_id=%s archived=%d deleted=%d delete_failed=%d",
                key.group_id(),
                len(entries),
                deleted,
                failed,
            )
        except Exception as exc:
            groups_map[key.group_id()] = _build_group_record(
                key=key,
                archive_path=archive_path,
                entries=entries,
                verified=False,
                verify_errors=[f"exception:{exc}"],
            )
            stats.groups_failed += 1
            stats.files_failed += len(entries)
            logger.exception(
                "[artifact-retention.group] status=failed group_id=%s reason=exception err=%s",
                key.group_id(),
                exc,
            )
            try:
                if tmp_tar.exists():
                    tmp_tar.unlink()
            except Exception:
                pass

    if mode == "minimal_disk" and not args.dry_run and stats.groups_failed == 0:
        extras = _collect_extra_minimal_disk_files(
            run_id=run_id,
            run_docked=run_docked,
            run_post_docked=run_post,
            pdb_id_filter=pdb_id_filter,
            include_globs=include_globs,
            exclude_globs=exclude_globs,
        )
        deleted, failed = _delete_files(extras)
        stats.extra_files_deleted += deleted
        stats.files_failed += failed
        logger.info(
            "[artifact-retention.minimal-disk] extras_deleted=%d extras_failed=%d",
            deleted,
            failed,
        )

    if not args.dry_run:
        _write_run_index(
            index_path,
            run_id=run_id,
            mode=mode,
            groups_map=groups_map,
        )

    logger.info("[artifact-retention.summary] %s", json.dumps(stats.as_dict(), sort_keys=True))
    return 0 if stats.groups_failed == 0 else 1


def _is_restore_target_safe(
    *,
    raw_target: Any,
    allowed_root_a: Path,
    allowed_root_b: Path,
) -> Tuple[bool, Optional[Path], str]:
    text = str(raw_target or "").strip()
    if not text:
        return False, None, "empty_path"
    target = Path(text).expanduser().resolve()
    if target.exists() and target.is_dir():
        return False, None, "target_is_directory"
    allowed_a = allowed_root_a.resolve()
    allowed_b = allowed_root_b.resolve()
    if target.is_relative_to(allowed_a) or target.is_relative_to(allowed_b):
        return True, target, "ok"
    return False, None, "outside_allowed_roots"


def _restore_archive(
    archive_path: Path,
    *,
    payload_label: str,
    payload: Mapping[str, Any],
    allowed_root_a: Path,
    allowed_root_b: Path,
    logger: logging.Logger,
) -> Tuple[int, int]:
    entries = payload.get("entries") or []
    if not isinstance(entries, list):
        raise RuntimeError(f"invalid restore entries: {payload_label}")

    restored = 0
    failed = 0
    with tempfile.TemporaryDirectory(prefix="artifact-retention-restore-") as tmp:
        tmp_root = Path(tmp)
        cmd = ["tar", "--zstd", "-xf", str(archive_path), "-C", str(tmp_root)]
        proc = subprocess.run(cmd, check=False)
        if proc.returncode != 0:
            raise RuntimeError(f"archive extraction failed for {archive_path} (code {proc.returncode})")

        for entry in entries:
            member_name = str(entry.get("member_name", ""))
            raw_original_path = entry.get("original_path", "")
            expected_sha = str(entry.get("sha256", ""))
            ok_target, original_path, reason = _is_restore_target_safe(
                raw_target=raw_original_path,
                allowed_root_a=allowed_root_a,
                allowed_root_b=allowed_root_b,
            )
            if not member_name or not ok_target or original_path is None:
                failed += 1
                logger.error(
                    "[artifact-retention.restore] status=failed member=%s reason=%s target=%s",
                    member_name or "(missing)",
                    reason,
                    raw_original_path,
                )
                continue
            src = tmp_root / member_name
            if not src.exists() or not src.is_file():
                failed += 1
                continue
            try:
                original_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, original_path)
                if expected_sha and _hash_file(original_path) != expected_sha:
                    failed += 1
                    logger.error(
                        "[artifact-retention.restore] status=failed member=%s reason=checksum_mismatch",
                        member_name,
                    )
                    continue
                restored += 1
            except Exception:
                failed += 1

    return restored, failed


def _resolve_archive_root_for_run(args: argparse.Namespace, run_id: str) -> Path:
    repo_root = (
        Path(args.repo_root).resolve()
        if args.repo_root
        else Path(__file__).resolve().parents[2]
    )
    cfg: Dict[str, Any] = {}
    try:
        cfg = dict(load_inputs())
    except Exception:
        cfg = {}
    _, _, run_docked, _ = _resolve_run_roots(
        args,
        repo_root=repo_root,
        cfg=cfg,
        run_id=run_id,
    )
    return run_docked / _ARCHIVE_ROOT_NAME


def _infer_run_id_from_archive_path(archive_path: Path) -> Optional[str]:
    parts = archive_path.resolve().parts
    for root_name in ("docked", "post_docked"):
        if root_name in parts:
            idx = parts.index(root_name)
            if idx + 1 < len(parts):
                token = str(parts[idx + 1]).strip()
                if token:
                    return token
    return None


def _resolve_restore_payload(
    args: argparse.Namespace,
    run_id: Optional[str],
) -> Tuple[Path, Dict[str, Any], str]:
    if args.restore_group:
        if not run_id:
            raise SystemExit("--restore-group requires --run-id")
        archive_root = _resolve_archive_root_for_run(args, run_id)
        index_path = archive_root / _INDEX_FILENAME
        if not index_path.exists():
            raise SystemExit(f"missing run-level index: {index_path}")
        index_payload = _load_run_index(index_path)
        for item in index_payload.get("groups", []):
            if str(item.get("group_id")) != str(args.restore_group):
                continue
            archive_path = _resolve_archive_path(item.get("archive_path", ""), archive_root)
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
                return archive_path, _load_manifest(manifest_path), str(manifest_path)
            raise SystemExit(f"group has no restore metadata: {args.restore_group}")
        raise SystemExit(f"group not found in index: {args.restore_group}")

    if not args.archive:
        raise SystemExit("--restore requires --archive <path> or --restore-group <group-id> with --run-id")

    archive_path = Path(args.archive).resolve()
    if args.manifest:
        manifest_path = Path(args.manifest).resolve()
        if not manifest_path.exists():
            raise SystemExit(f"manifest not found: {manifest_path}")
        return archive_path, _load_manifest(manifest_path), str(manifest_path)

    candidate_run_id = run_id or _infer_run_id_from_archive_path(archive_path)
    if candidate_run_id:
        archive_root = _resolve_archive_root_for_run(args, candidate_run_id)
        index_path = archive_root / _INDEX_FILENAME
        if index_path.exists():
            index_payload = _load_run_index(index_path)
            archive_norm = archive_path.resolve()
            for item in index_payload.get("groups", []):
                item_archive = _resolve_archive_path(item.get("archive_path", ""), archive_root)
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
                        return archive_path, _load_manifest(manifest_path), str(manifest_path)

    legacy_manifest = archive_path.with_name(_MANIFEST_FILENAME)
    if legacy_manifest.exists():
        return archive_path, _load_manifest(legacy_manifest), str(legacy_manifest)
    raise SystemExit(
        f"unable to resolve restore metadata for archive: {archive_path}; "
        "use --manifest or --restore-group with --run-id"
    )


def _filter_restore_entries(
    entries: Sequence[Dict[str, Any]],
    *,
    stage: Optional[str],
    mode: Optional[str],
    path_globs: Sequence[str],
) -> List[Dict[str, Any]]:
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
            matched = False
            for pat in path_globs:
                if fnmatch.fnmatch(member, pat) or fnmatch.fnmatch(original, pat):
                    matched = True
                    break
            if not matched:
                continue
        out.append(entry)
    return out


def _run_restore(args: argparse.Namespace) -> int:
    logger = logging.getLogger("artifact-retention")
    preflight_ok, preflight_reason = _tooling_preflight()
    if not preflight_ok:
        logger.error(
            "[artifact-retention.restore] reason=tooling_preflight_failed detail=%s",
            preflight_reason,
        )
        return 1

    run_id = str(args.run_id or "").strip() or None
    archive_path, payload, payload_label = _resolve_restore_payload(args, run_id=run_id)
    if not archive_path.exists():
        raise SystemExit(f"archive not found: {archive_path}")
    payload_run_id = _normalize_run_id(payload.get("group", {}).get("run_id") or "")
    if not payload_run_id:
        payload_run_id = _normalize_run_id(
            (payload.get("entries") or [{}])[0].get("run_id", "")
        )
    effective_run_id = run_id or payload_run_id
    if not effective_run_id:
        raise SystemExit("manifest missing run_id; cannot validate restore target roots")

    repo_root = (
        Path(args.repo_root).resolve()
        if args.repo_root
        else Path(__file__).resolve().parents[2]
    )
    cfg: Dict[str, Any] = {}
    try:
        cfg = dict(load_inputs())
    except Exception:
        cfg = {}
    _, _, run_docked, run_post = _resolve_run_roots(
        args,
        repo_root=repo_root,
        cfg=cfg,
        run_id=effective_run_id,
    )

    all_entries = payload.get("entries") or []
    if not isinstance(all_entries, list):
        raise SystemExit("restore payload has invalid entries")
    selected_entries = _filter_restore_entries(
        all_entries,
        stage=args.restore_stage,
        mode=args.restore_mode,
        path_globs=list(args.restore_path_glob or []),
    )
    if not selected_entries:
        raise SystemExit("no entries matched restore filters")
    restore_payload = dict(payload)
    restore_payload["entries"] = selected_entries

    restored, failed = _restore_archive(
        archive_path,
        payload_label=payload_label,
        payload=restore_payload,
        allowed_root_a=run_docked,
        allowed_root_b=run_post,
        logger=logger,
    )
    payload = {
        "action": "restore",
        "archive": str(archive_path),
        "restore_source": payload_label,
        "entries_selected": len(selected_entries),
        "restored": restored,
        "failed": failed,
    }
    logger.info("[artifact-retention.summary] %s", json.dumps(payload, sort_keys=True))
    return 0 if failed == 0 else 1


def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Archive bulky per-ligand artifacts under docked/post_docked into deterministic "
            ".tar.zst groups with verification and restore support."
        )
    )
    ap.add_argument("--run-id", help="Run ID under docked/ and post_docked/")
    ap.add_argument(
        "--mode",
        choices=["rerun_safe", "minimal_disk", "off"],
        default="rerun_safe",
        help="Retention mode (default: rerun_safe)",
    )
    ap.add_argument("--pdb-id", help="Optional PDB scope filter (e.g., 1ABC)")
    ap.add_argument("--dry-run", action="store_true", help="Plan actions only; no archive/deletion")
    ap.add_argument("--overwrite", action="store_true", help="Rebuild archives even if already verified")
    ap.add_argument(
        "--include-glob",
        action="append",
        default=[],
        help="Glob include filter matched against archive member path (repeatable)",
    )
    ap.add_argument(
        "--exclude-glob",
        action="append",
        default=[],
        help="Glob exclude filter matched against archive member path (repeatable)",
    )
    ap.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1), help="Worker count for hashing/compression")
    ap.add_argument("--verify-only", action="store_true", help="Verify existing archives from run index/manifests only")
    ap.add_argument("--restore", action="store_true", help="Restore files from one archive group")
    ap.add_argument("--archive", help="Archive file path for --restore")
    ap.add_argument("--manifest", help="Manifest JSON path for --restore (optional with --archive)")
    ap.add_argument("--restore-group", help="Group id from archive_index.json for --restore")
    ap.add_argument("--restore-stage", help="Optional stage_dir filter when restoring")
    ap.add_argument("--restore-mode", help="Optional mode filter when restoring")
    ap.add_argument(
        "--restore-path-glob",
        action="append",
        default=[],
        help="Optional glob filter for member_name/original_path when restoring (repeatable)",
    )
    ap.add_argument("--repo-root", help="Override repository root (defaults from this module path)")
    ap.add_argument("--docked-root", help="Override docked root (run-id appended unless already present)")
    ap.add_argument("--post-docked-root", help="Override post_docked root (run-id appended unless already present)")
    ap.add_argument("--verbose", action="store_true", help="Verbose logging")
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = _build_parser()
    args = ap.parse_args(argv)
    _setup_logging(verbose=bool(args.verbose))

    if int(args.threads) <= 0:
        raise SystemExit("--threads must be > 0")

    if args.restore:
        return _run_restore(args)
    return _run_retention(args)


if __name__ == "__main__":
    raise SystemExit(main())
