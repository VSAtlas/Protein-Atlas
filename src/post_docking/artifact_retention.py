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
    stage_dir: str
    mode: str

    def as_dict(self) -> Dict[str, str]:
        return {
            "run_id": self.run_id,
            "source_root": self.source_root,
            "pdb_id": self.pdb_id,
            "variant": self.variant,
            "ph": self.ph,
            "stage_dir": self.stage_dir,
            "mode": self.mode,
        }

    def group_id(self) -> str:
        return (
            f"run={self.run_id}|source={self.source_root}|pdb={self.pdb_id}|"
            f"variant={self.variant}|ph={self.ph}|stage={self.stage_dir}|mode={self.mode}"
        )


@dataclass(frozen=True)
class CandidateFile:
    key: GroupKey
    original_path: Path
    member_name: str
    file_type: str


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

    stage_dir = "root"
    for piece in parts[idx:-1]:
        if piece.lower().startswith("stage"):
            stage_dir = piece
            break
    if stage_dir == "root" and len(parts) >= 2:
        stage_dir = parts[-2]

    mode = _infer_mode(parts)
    return GroupKey(
        run_id=run_id,
        source_root=source_root,
        pdb_id=pdb_id,
        variant=variant,
        ph=ph,
        stage_dir=stage_dir,
        mode=mode,
    )


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
    if threads <= 1 or len(paths) <= 1:
        return {p: _hash_file(p) for p in paths}

    out: Dict[Path, str] = {}
    with ThreadPoolExecutor(max_workers=threads) as pool:
        for p, sha in zip(paths, pool.map(_hash_file, paths)):
            out[p] = sha
    return out


def _archive_group_paths(archive_root: Path, key: GroupKey) -> Tuple[Path, Path]:
    group_dir = (
        archive_root
        / _safe_token(key.source_root)
        / _safe_token(key.pdb_id)
        / _safe_token(key.variant)
        / _safe_token(key.ph)
        / _safe_token(key.stage_dir)
        / _safe_token(key.mode)
    )
    return group_dir / _ARCHIVE_FILENAME, group_dir / _MANIFEST_FILENAME


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


def _write_manifest(
    manifest_path: Path,
    *,
    key: GroupKey,
    archive_path: Path,
    entries: Sequence[Dict[str, Any]],
    verified: bool,
    verify_errors: Sequence[str],
) -> None:
    payload = {
        "schema_version": 1,
        "generated_at": _utc_now(),
        "group": key.as_dict(),
        "group_id": key.group_id(),
        "archive_path": str(archive_path),
        "manifest_path": str(manifest_path),
        "entry_count": len(entries),
        "verified": bool(verified),
        "verify_errors": list(verify_errors),
        "entries": list(entries),
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=False), encoding="utf-8")


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


def _refresh_run_index(archive_root: Path, run_id: str, mode: str) -> None:
    manifests = sorted(archive_root.rglob(_MANIFEST_FILENAME))
    groups: List[Dict[str, Any]] = []
    for manifest_path in manifests:
        try:
            data = _load_manifest(manifest_path)
        except Exception:
            continue
        group = data.get("group") or {}
        groups.append(
            {
                "group_id": data.get("group_id"),
                "group": group,
                "archive_path": data.get("archive_path"),
                "manifest_path": str(manifest_path),
                "entry_count": data.get("entry_count", 0),
                "verified": bool(data.get("verified", False)),
            }
        )

    payload = {
        "schema_version": 1,
        "generated_at": _utc_now(),
        "run_id": run_id,
        "mode": mode,
        "group_count": len(groups),
        "groups": groups,
    }
    (archive_root / _INDEX_FILENAME).write_text(
        json.dumps(payload, indent=2, sort_keys=False),
        encoding="utf-8",
    )


def _build_entries(
    files: Sequence[CandidateFile],
    archive_path: Path,
    threads: int,
) -> List[Dict[str, Any]]:
    paths = [f.original_path for f in files]
    hashes = _hash_files(paths, threads=threads)
    entries: List[Dict[str, Any]] = []
    for f in files:
        st = f.original_path.stat()
        entry = {
            "run_id": f.key.run_id,
            "pdb_id": f.key.pdb_id,
            "variant": f.key.variant,
            "ph": f.key.ph,
            "stage_dir": f.key.stage_dir,
            "mode": f.key.mode,
            "original_path": str(f.original_path),
            "archive_path": str(archive_path),
            "member_name": f.member_name,
            "size_bytes": int(st.st_size),
            "mtime": float(st.st_mtime),
            "sha256": hashes[f.original_path],
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

    if args.verify_only:
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
        archive_path, manifest_path = _archive_group_paths(archive_root, key)
        stats.groups_planned += 1
        discovered_members = {item.member_name for item in files}

        if (
            not args.overwrite
            and archive_path.exists()
            and manifest_path.exists()
        ):
            try:
                payload = _load_manifest(manifest_path)
                manifest_entries = payload.get("entries") or []
                manifest_members = {
                    str(entry.get("member_name", ""))
                    for entry in manifest_entries
                    if isinstance(entry, dict)
                }
                manifest_matches = (
                    len(manifest_entries) == len(files)
                    and manifest_members == discovered_members
                )
                if bool(payload.get("verified", False)) and manifest_matches:
                    stats.groups_skipped += 1
                    stats.files_skipped += len(files)
                    logger.info(
                        "[artifact-retention.group] action=skip group_id=%s reason=already_archived_verified files=%d",
                        key.group_id(),
                        len(files),
                    )
                    continue
                if bool(payload.get("verified", False)) and not manifest_matches:
                    logger.warning(
                        "[artifact-retention.group] action=rebuild group_id=%s reason=manifest_mismatch discovered=%d manifest=%d",
                        key.group_id(),
                        len(files),
                        len(manifest_entries),
                    )
            except Exception:
                logger.warning(
                    "[artifact-retention.group] action=rebuild group_id=%s reason=manifest_unreadable",
                    key.group_id(),
                )

        entries = _build_entries(files, archive_path=archive_path, threads=max(1, int(args.threads)))
        path_by_member = {item.member_name: item.original_path for item in files}

        if args.dry_run:
            logger.info(
                "[artifact-retention.group] action=dry_run_plan group_id=%s files=%d archive=%s manifest=%s",
                key.group_id(),
                len(files),
                archive_path,
                manifest_path,
            )
            continue

        archive_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_tar = archive_path.with_suffix(".tar.tmp")
        if tmp_tar.exists():
            tmp_tar.unlink()
        if archive_path.exists() and args.overwrite:
            archive_path.unlink()
        if manifest_path.exists() and args.overwrite:
            manifest_path.unlink()

        try:
            _write_tar(entries=entries, path_by_member=path_by_member, tar_path=tmp_tar)
            _compress_tar_to_zst(tmp_tar, archive_path=archive_path, threads=max(1, int(args.threads)))
            ok, verify_errors = _verify_archive(archive_path, entries)
            _write_manifest(
                manifest_path,
                key=key,
                archive_path=archive_path,
                entries=entries,
                verified=ok,
                verify_errors=verify_errors,
            )
            if not ok:
                stats.groups_failed += 1
                stats.files_failed += len(files)
                logger.error(
                    "[artifact-retention.group] status=failed group_id=%s reason=verify_failed errors=%s",
                    key.group_id(),
                    ";".join(verify_errors[:10]),
                )
                continue

            stats.groups_archived += 1
            stats.groups_verified += 1
            stats.files_archived += len(files)
            to_delete = [item.original_path for item in files]
            deleted, failed = _delete_files(to_delete)
            stats.files_deleted += deleted
            stats.files_failed += failed
            logger.info(
                "[artifact-retention.group] status=ok group_id=%s archived=%d deleted=%d delete_failed=%d",
                key.group_id(),
                len(files),
                deleted,
                failed,
            )
        except Exception as exc:
            stats.groups_failed += 1
            stats.files_failed += len(files)
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
        _refresh_run_index(archive_root, run_id=run_id, mode=mode)

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
    manifest_path: Path,
    *,
    payload: Mapping[str, Any],
    allowed_root_a: Path,
    allowed_root_b: Path,
    logger: logging.Logger,
) -> Tuple[int, int]:
    entries = payload.get("entries") or []
    if not isinstance(entries, list):
        raise RuntimeError(f"invalid manifest entries: {manifest_path}")

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


def _resolve_restore_target(args: argparse.Namespace, run_id: Optional[str]) -> Tuple[Path, Path]:
    if args.archive:
        archive_path = Path(args.archive).resolve()
        manifest_path = archive_path.with_suffix("").with_suffix(".manifest.json")
        if archive_path.name.endswith(".tar.zst"):
            manifest_path = archive_path.with_name(_MANIFEST_FILENAME)
        if args.manifest:
            manifest_path = Path(args.manifest).resolve()
        return archive_path, manifest_path

    if not args.restore_group or not run_id:
        raise SystemExit("--restore requires --archive <path> or --restore-group <group-id> with --run-id")

    repo_root = Path(args.repo_root).resolve() if args.repo_root else Path(__file__).resolve().parents[2]
    cfg: Dict[str, Any] = {}
    try:
        cfg = dict(load_inputs())
    except Exception:
        cfg = {}

    overall = Path(str(cfg.get("OVERALL_DIR", repo_root))).resolve()
    docked_base = (
        Path(args.docked_root).resolve()
        if args.docked_root
        else Path(str(cfg.get("DOCKED_DIR", overall / "docked"))).resolve()
    )
    run_docked = docked_base if docked_base.name == run_id else docked_base / run_id
    archive_root = run_docked / _ARCHIVE_ROOT_NAME
    index_path = archive_root / _INDEX_FILENAME
    if not index_path.exists():
        raise SystemExit(f"missing run-level index: {index_path}")
    index_payload = json.loads(index_path.read_text(encoding="utf-8"))
    for item in index_payload.get("groups", []):
        if str(item.get("group_id")) == str(args.restore_group):
            archive_path = Path(str(item.get("archive_path"))).resolve()
            manifest_path = Path(str(item.get("manifest_path"))).resolve()
            return archive_path, manifest_path
    raise SystemExit(f"group not found in index: {args.restore_group}")


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
    archive_path, manifest_path = _resolve_restore_target(args, run_id=run_id)
    if not archive_path.exists():
        raise SystemExit(f"archive not found: {archive_path}")
    if not manifest_path.exists():
        raise SystemExit(f"manifest not found: {manifest_path}")
    payload = _load_manifest(manifest_path)
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

    restored, failed = _restore_archive(
        archive_path,
        manifest_path,
        payload=payload,
        allowed_root_a=run_docked,
        allowed_root_b=run_post,
        logger=logger,
    )
    payload = {
        "action": "restore",
        "archive": str(archive_path),
        "manifest": str(manifest_path),
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
    ap.add_argument("--verify-only", action="store_true", help="Verify existing archives/manifests only")
    ap.add_argument("--restore", action="store_true", help="Restore files from one archive group")
    ap.add_argument("--archive", help="Archive file path for --restore")
    ap.add_argument("--manifest", help="Manifest JSON path for --restore (optional with --archive)")
    ap.add_argument("--restore-group", help="Group id from archive_index.json for --restore")
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
