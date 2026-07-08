# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import json
import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from post_docking.artifact_retention_support import (
    _MANIFEST_FILENAME,
    _ARCHIVE_ROOT_NAME,
    _INDEX_FILENAME,
    GroupKey,
    RetentionStats,
    archive_group_path as _archive_group_path,
    build_entries as _build_entries,
    build_group_record as _build_group_record,
    collect_extra_minimal_disk_files as _collect_extra_minimal_disk_files,
    compress_tar_to_zst as _compress_tar_to_zst,
    delete_files as _delete_files,
    discover_candidates as _discover_candidates,
    filter_restore_entries as _filter_restore_entries,
    load_inputs_safe as _load_inputs_safe,
    load_manifest as _load_manifest,
    load_run_index as _load_run_index,
    normalize_pdb_id as _normalize_pdb_id,
    normalize_ph_token as _normalize_ph_token,
    normalize_ph_filter as _normalize_ph_filter,
    normalize_retention_mode as _normalize_retention_mode,
    normalize_run_id as _normalize_run_id,
    normalize_variant_filter as _normalize_variant_filter,
    normalize_variant_token as _normalize_variant_token,
    resolve_archive_path as _resolve_archive_path,
    resolve_combo_scope as _resolve_combo_scope,
    resolve_restore_payload as _resolve_restore_payload,
    resolve_run_roots as _resolve_run_roots,
    restore_archive as _restore_archive,
    setup_logging as _setup_logging,
    tooling_preflight as _tooling_preflight,
    verify_archive as _verify_archive,
    write_run_index as _write_run_index,
    write_tar as _write_tar,
)

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
    variant_filter = _normalize_variant_filter(args.variant)
    ph_filter = _normalize_ph_filter(args.ph)
    combo_scope = _resolve_combo_scope(args)
    if combo_scope and (pdb_id_filter or variant_filter or ph_filter):
        raise SystemExit("--combo/--combo-file is mutually exclusive with --pdb-id/--variant/--ph")
    if (variant_filter or ph_filter) and not pdb_id_filter:
        raise SystemExit("--variant/--ph filters require --pdb-id")
    try:
        group_workers = max(1, int(args.group_workers))
    except Exception:
        raise SystemExit("--group-workers must be > 0")
    threads_total = max(1, int(args.threads))
    threads_per_group = max(1, threads_total // group_workers)

    repo_root = (
        Path(args.repo_root).resolve()
        if args.repo_root
        else Path(__file__).resolve().parents[2]
    )
    cfg = _load_inputs_safe(logger=logger)
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
        "[artifact-retention.start] run_id=%s pdb_id=%s variant=%s ph=%s combo_count=%d mode=%s dry_run=%s verify_only=%s overwrite=%s docked_root=%s post_docked_root=%s",
        run_id,
        pdb_id_filter or "*",
        variant_filter or "*",
        ph_filter or "*",
        len(combo_scope),
        mode,
        str(bool(args.dry_run)).lower(),
        str(bool(args.verify_only)).lower(),
        str(bool(args.overwrite)).lower(),
        run_docked,
        run_post,
    )
    logger.info(
        "[artifact-retention.batch] run_id=%s combo_count=%d group_workers=%d threads_total=%d threads_per_group=%d",
        run_id,
        len(combo_scope),
        group_workers,
        threads_total,
        threads_per_group,
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
        variant_filter=variant_filter,
        ph_filter=ph_filter,
        include_globs=include_globs,
        exclude_globs=exclude_globs,
        logger=logger,
    )
    if combo_scope:
        groups = {
            key: files
            for key, files in groups.items()
            if (
                str(key.pdb_id).upper(),
                _normalize_variant_token(str(key.variant)),
                _normalize_ph_token(str(key.ph)),
            )
            in combo_scope
        }
    stats.groups_discovered = len(groups)
    stats.files_discovered = sum(len(v) for v in groups.values())

    sorted_group_keys = sorted(groups.keys(), key=lambda k: k.group_id())
    stats.groups_planned = len(sorted_group_keys)
    retention_wall_start = datetime.now(timezone.utc)

    def _process_group(key: GroupKey) -> Dict[str, Any]:
        files = groups[key]
        archive_path = _archive_group_path(archive_root, key)
        result: Dict[str, Any] = {
            "group_id": key.group_id(),
            "record": None,
            "groups_skipped": 0,
            "groups_archived": 0,
            "groups_verified": 0,
            "groups_failed": 0,
            "files_skipped": 0,
            "files_archived": 0,
            "files_deleted": 0,
            "files_failed": 0,
        }
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
                result["groups_skipped"] = 1
                result["files_skipped"] = len(files)
                logger.info(
                    "[artifact-retention.group] action=skip group_id=%s reason=already_archived_verified files=%d",
                    key.group_id(),
                    len(files),
                )
                return result
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
            return result

        entries = _build_entries(
            files,
            archive_path=archive_path,
            threads=threads_per_group,
        )
        if len(entries) != len(files):
            vanished = len(files) - len(entries)
            if vanished > 0:
                result["files_skipped"] += vanished
                logger.warning(
                    "[artifact-retention.group] action=partial group_id=%s reason=files_missing_during_hash_or_stat discovered=%d stable=%d",
                    key.group_id(),
                    len(files),
                    len(entries),
                )
        if not entries:
            result["groups_skipped"] = 1
            logger.info(
                "[artifact-retention.group] action=skip group_id=%s reason=no_stable_files",
                key.group_id(),
            )
            return result

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
            _compress_tar_to_zst(
                tmp_tar,
                archive_path=archive_path,
                threads=threads_per_group,
            )
            ok, verify_errors = _verify_archive(archive_path, entries)
            result["record"] = _build_group_record(
                key=key,
                archive_path=archive_path,
                entries=entries,
                verified=ok,
                verify_errors=verify_errors,
            )
            if not ok:
                result["groups_failed"] = 1
                result["files_failed"] = len(entries)
                logger.error(
                    "[artifact-retention.group] status=failed group_id=%s reason=verify_failed errors=%s",
                    key.group_id(),
                    ";".join(verify_errors[:10]),
                )
                return result

            result["groups_archived"] = 1
            result["groups_verified"] = 1
            result["files_archived"] = len(entries)
            deleted, failed = _delete_files(list(path_by_member.values()))
            result["files_deleted"] = deleted
            result["files_failed"] = failed
            logger.info(
                "[artifact-retention.group] status=ok group_id=%s archived=%d deleted=%d delete_failed=%d",
                key.group_id(),
                len(entries),
                deleted,
                failed,
            )
            return result
        except Exception as exc:
            result["record"] = _build_group_record(
                key=key,
                archive_path=archive_path,
                entries=entries,
                verified=False,
                verify_errors=[f"exception:{exc}"],
            )
            result["groups_failed"] = 1
            result["files_failed"] = len(entries)
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
            return result

    completed = 0
    total = len(sorted_group_keys)
    if group_workers == 1 or total <= 1:
        for key in sorted_group_keys:
            payload = _process_group(key)
            record = payload.get("record")
            if isinstance(record, dict):
                groups_map[str(payload.get("group_id", key.group_id()))] = record
            for stat_key in (
                "groups_skipped",
                "groups_archived",
                "groups_verified",
                "groups_failed",
                "files_skipped",
                "files_archived",
                "files_deleted",
                "files_failed",
            ):
                stats.__dict__[stat_key] += int(payload.get(stat_key, 0))
            completed += 1
            logger.info(
                "[artifact-retention.progress] completed=%d total=%d failed=%d",
                completed,
                total,
                stats.groups_failed,
            )
    else:
        with ThreadPoolExecutor(max_workers=group_workers) as pool:
            future_map = {pool.submit(_process_group, key): key for key in sorted_group_keys}
            for fut in as_completed(future_map):
                key = future_map[fut]
                try:
                    payload = fut.result()
                except Exception as exc:
                    payload = {
                        "group_id": key.group_id(),
                        "record": None,
                        "groups_failed": 1,
                        "files_failed": len(groups.get(key, [])),
                    }
                    logger.exception(
                        "[artifact-retention.group] status=failed group_id=%s reason=future_exception err=%s",
                        key.group_id(),
                        exc,
                    )
                record = payload.get("record")
                if isinstance(record, dict):
                    groups_map[str(payload.get("group_id", key.group_id()))] = record
                for stat_key in (
                    "groups_skipped",
                    "groups_archived",
                    "groups_verified",
                    "groups_failed",
                    "files_skipped",
                    "files_archived",
                    "files_deleted",
                    "files_failed",
                ):
                    stats.__dict__[stat_key] += int(payload.get(stat_key, 0))
                completed += 1
                logger.info(
                    "[artifact-retention.progress] completed=%d total=%d failed=%d",
                    completed,
                    total,
                    stats.groups_failed,
                )

    elapsed_sec = max(
        0.0,
        (datetime.now(timezone.utc) - retention_wall_start).total_seconds(),
    )
    logger.info(
        "[artifact-retention.batch.done] run_id=%s groups_total=%d groups_failed=%d elapsed_sec=%.2f",
        run_id,
        total,
        stats.groups_failed,
        elapsed_sec,
    )

    if mode == "minimal_disk" and not args.dry_run and stats.groups_failed == 0:
        extras = _collect_extra_minimal_disk_files(
            run_id=run_id,
            run_docked=run_docked,
            run_post_docked=run_post,
            pdb_id_filter=pdb_id_filter,
            variant_filter=variant_filter,
            ph_filter=ph_filter,
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
    cfg = _load_inputs_safe(logger=logger)
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
    ap.add_argument("--variant", help="Optional variant scope filter (e.g., HOLO/APO)")
    ap.add_argument("--ph", help="Optional pH scope filter (e.g., pH7_7)")
    ap.add_argument(
        "--combo",
        action="append",
        default=[],
        help="Optional combo scope token PDB:VARIANT:PH (repeatable)",
    )
    ap.add_argument(
        "--combo-file",
        help="Optional newline-delimited combo scope file (PDB:VARIANT:PH per line)",
    )
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
    ap.add_argument(
        "--group-workers",
        type=int,
        default=1,
        help="Concurrent group workers (default: 1)",
    )
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
    if int(args.group_workers) <= 0:
        raise SystemExit("--group-workers must be > 0")

    if args.restore:
        return _run_restore(args)
    return _run_retention(args)


if __name__ == "__main__":
    raise SystemExit(main())
