from __future__ import annotations

import argparse
import base64
import csv
import html
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from cli.postrun_hooks_support import _repo_root
from config.output_paths import run_output_dir
from config.runtime_config import load_config
from docking.capture_pose import capture_pose


@dataclass(frozen=True)
class ScreenshotSelection:
    pdb_id: str
    ligand: str
    score: float | None
    score_source: str
    pose_path: Path
    receptor_path: Path
    output_prefix: Path
    selection_label: str = ""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="atlas screenshot",
        description="Render active-site screenshots for selected docked ligand poses.",
    )
    parser.add_argument("run_id", help="Atlas run ID under outputs/docked/<RUN_ID>/")
    parser.add_argument("--pdb", action="append", default=[], help="PDB ID to render")
    parser.add_argument("--ligand", action="append", default=[], help="Ligand name/stem selector")
    parser.add_argument("--top", type=int, default=0, help="Render top N ligands per PDB")
    parser.add_argument(
        "--representatives",
        action="store_true",
        help="Render representative best, median, and worst scored ligands per PDB.",
    )
    parser.add_argument(
        "--gallery-html",
        action="store_true",
        help="Write an HTML review gallery next to generated screenshots.",
    )
    parser.add_argument("--stage", default="stage3", help="Docking stage directory to render")
    parser.add_argument("--variant", help="Optional variant filter for archived outputs")
    parser.add_argument("--ph", help="Optional pH label filter for archived outputs")
    parser.add_argument(
        "--without-ligand",
        action="store_true",
        help="Use the ligand to define the active site but hide it in the screenshots.",
    )
    output_group = parser.add_mutually_exclusive_group()
    output_group.add_argument(
        "--out-dir",
        help="Output directory; defaults to outputs/data/<RUN_ID>/screenshots/",
    )
    output_group.add_argument(
        "--output-subdir",
        help=(
            "Safe relative subdirectory beneath "
            "outputs/data/<RUN_ID>/screenshots/"
        ),
    )
    parser.add_argument("--label-top-n-res", type=int, default=0)
    parser.add_argument("--cutoff", type=float, default=6.0)
    parser.add_argument(
        "--protein-transparency",
        type=float,
        default=0.65,
        help="Surface transparency; 1.0 is invisible and 0.0 is opaque.",
    )
    parser.add_argument(
        "--protein-style",
        choices=("cartoon", "lines", "mesh", "surface", "both", "none"),
        default="cartoon",
        help="Protein pocket representation for ligand inspection.",
    )
    parser.add_argument(
        "--surface-carve-cutoff",
        type=float,
        default=5.0,
        help="PyMOL surface-carve distance around the ligand.",
    )
    parser.add_argument(
        "--surface-carve-normal-cutoff",
        type=float,
        default=-0.1,
        help="Hide carved surface triangles not facing the ligand.",
    )
    parser.add_argument(
        "--clip-slab",
        type=float,
        default=0.0,
        help="Camera slab thickness around the ligand-pocket context; 0 disables.",
    )
    parser.add_argument(
        "--view-context",
        choices=("full", "pocket"),
        default="full",
        help="Camera framing context; full shows the protein, pocket zooms to the ligand site.",
    )
    parser.add_argument(
        "--zoom-buffer",
        type=float,
        default=4.0,
        help="Angstrom buffer around the rendered pocket context.",
    )
    parser.add_argument("--viewport", default="1600x1200")
    parser.add_argument("--ray", action="store_true", help="Ray-trace PNGs for publication output.")
    parser.add_argument(
        "--no-restore",
        action="store_true",
        help="Do not restore archived poses; only render files already on disk.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print selected receptor/pose/output paths without rendering.",
    )
    return parser


def run_screenshot(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(list(argv or []))
    _validate_args(args)
    selections = _select_screenshots(args)
    if not selections:
        print("atlas screenshot: no matching poses found", file=sys.stderr)
        return 1
    if args.dry_run:
        _print_selections(selections)
        return 0
    return _render_selections(selections, args)


def _validate_args(args: argparse.Namespace) -> None:
    if not args.pdb:
        raise SystemExit("atlas screenshot requires at least one --pdb")
    if not args.ligand and int(args.top or 0) <= 0 and not bool(args.representatives):
        raise SystemExit("provide --ligand <name>, --top N, or --representatives")


def _select_screenshots(args: argparse.Namespace) -> list[ScreenshotSelection]:
    repo_root = _repo_root()
    cfg = load_config("config.txt", base_dir=repo_root)
    out_dir = _resolve_output_dir(repo_root, args)
    return _build_selections(repo_root, cfg, args, out_dir)


def _render_selections(
    selections: Sequence[ScreenshotSelection],
    args: argparse.Namespace,
) -> int:
    out_dir = _resolve_output_dir(_repo_root(), args)
    viewport = _parse_viewport(args.viewport)
    rendered = 0
    failed = 0
    rendered_selections: list[ScreenshotSelection] = []
    for selection in selections:
        selection.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        try:
            capture_pose(
                selection.receptor_path,
                selection.pose_path,
                selection.output_prefix,
                top_n_residues=int(args.label_top_n_res),
                proximity_cutoff=float(args.cutoff),
                viewport=viewport,
                zoom_buffer=float(args.zoom_buffer),
                protein_transparency=float(args.protein_transparency),
                protein_style=str(args.protein_style),
                view_context=str(args.view_context),
                surface_carve_cutoff=float(args.surface_carve_cutoff),
                surface_carve_normal_cutoff=float(args.surface_carve_normal_cutoff),
                clip_slab=float(args.clip_slab),
                ray=bool(args.ray),
                hide_ligand=bool(args.without_ligand),
            )
            rendered += 1
            rendered_selections.append(selection)
            print(f"rendered: {selection.output_prefix}_front.png")
        except Exception as exc:
            failed += 1
            print(
                f"atlas screenshot: render failed ligand={selection.ligand} err={exc}",
                file=sys.stderr,
            )
    print(f"screenshots: rendered={rendered} failed={failed} out={out_dir}")
    if args.gallery_html and rendered_selections:
        gallery = _write_gallery_html(rendered_selections, args, out_dir)
        print(f"gallery: {gallery}")
    return 0 if failed == 0 else 1


def _resolve_output_dir(repo_root: Path, args: argparse.Namespace) -> Path:
    if args.out_dir:
        return Path(args.out_dir).expanduser().resolve()
    screenshots_root = run_output_dir(repo_root, "data", args.run_id) / "screenshots"
    output_subdir = str(getattr(args, "output_subdir", "") or "").strip()
    if not output_subdir:
        return screenshots_root
    relative = Path(output_subdir)
    if relative.is_absolute() or ".." in relative.parts:
        raise SystemExit("--output-subdir must be relative and must not contain '..'")
    root_resolved = screenshots_root.resolve()
    candidate = (root_resolved / relative).resolve()
    try:
        candidate.relative_to(root_resolved)
    except ValueError as exc:
        raise SystemExit(
            "--output-subdir must remain beneath the run screenshots directory"
        ) from exc
    return candidate


def _build_selections(
    repo_root: Path,
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
    out_dir: Path,
) -> list[ScreenshotSelection]:
    selections: list[ScreenshotSelection] = []
    for pdb_id in _clean_tokens(args.pdb):
        receptor = _resolve_receptor(repo_root, cfg, args.run_id, pdb_id)
        if receptor is None:
            print(f"atlas screenshot: receptor not found for {pdb_id}", file=sys.stderr)
            continue
        candidates = _candidate_ligands(repo_root, args, pdb_id)
        if not candidates:
            print(f"atlas screenshot: no score rows for {pdb_id}", file=sys.stderr)
            continue
        for ligand, score, source, label in candidates:
            pose = _pose_path(repo_root, args.run_id, pdb_id, ligand, args.stage)
            if not pose.exists() and not args.no_restore and not args.dry_run:
                _restore_pose_from_archive(repo_root, args, pose)
            if not pose.exists():
                print(f"atlas screenshot: pose not found: {pose}", file=sys.stderr)
                continue
            selections.append(
                _make_selection(
                    pdb_id=pdb_id,
                    ligand=ligand,
                    score=score,
                    score_source=source,
                    pose_path=pose,
                    receptor_path=receptor,
                    out_dir=out_dir,
                    without_ligand=bool(args.without_ligand),
                    selection_label=label,
                )
            )
    return selections


def _clean_tokens(values: Iterable[str]) -> list[str]:
    out: list[str] = []
    for value in values:
        token = str(value or "").strip()
        if token:
            out.append(token.upper())
    return out


def _resolve_receptor(repo_root: Path, cfg: Mapping[str, Any], run_id: str, pdb_id: str) -> Path | None:
    output_root = Path(str(cfg.get("OUTPUT_DIR") or run_output_dir(repo_root, "processed_pdbs", "")))
    input_root = Path(str(cfg.get("INPUT_DIR") or repo_root / "input_pdbs"))
    run_root = run_output_dir(repo_root, "processed_pdbs", run_id)
    candidates = [
        run_root / pdb_id / "receptor" / f"{pdb_id}_cleaned.pdb",
        run_root / pdb_id / "receptor" / f"{pdb_id}.pdb",
        run_root / pdb_id / "HOLO" / "receptor" / f"{pdb_id}_cleaned.pdb",
        run_root / pdb_id / "APO" / "receptor" / f"{pdb_id}_cleaned.pdb",
        output_root / pdb_id / "receptor" / f"{pdb_id}_cleaned.pdb",
        output_root / pdb_id / "receptor" / f"{pdb_id}.pdb",
        output_root / pdb_id / "HOLO" / "receptor" / f"{pdb_id}_cleaned.pdb",
        output_root / pdb_id / "APO" / "receptor" / f"{pdb_id}_cleaned.pdb",
        run_output_dir(repo_root, "processed_pdbs", pdb_id) / "receptor" / f"{pdb_id}_cleaned.pdb",
        input_root / f"{pdb_id}.pdb",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return None


def _candidate_ligands(
    repo_root: Path,
    args: argparse.Namespace,
    pdb_id: str,
) -> list[tuple[str, float | None, str, str]]:
    explicit = [str(item).strip() for item in args.ligand if str(item).strip()]
    if explicit:
        return [(item, None, "explicit", "explicit") for item in explicit]
    top_n = max(0, int(args.top or 0))
    scored = _read_consensus_ligands(repo_root, args.run_id, pdb_id)
    if not scored:
        scored = _read_docking_ligands(repo_root, args.run_id, pdb_id, args.stage)
    if args.representatives:
        return _representative_ligands(scored)
    if top_n <= 0:
        return []
    return [(ligand, score, source, f"top_{idx}") for idx, (ligand, score, source) in enumerate(scored[:top_n], 1)]


def _representative_ligands(
    scored: Sequence[tuple[str, float | None, str]],
) -> list[tuple[str, float | None, str, str]]:
    if not scored:
        return []
    indexes = [0, len(scored) // 2, len(scored) - 1]
    labels = ["best", "median", "worst"]
    out: list[tuple[str, float | None, str, str]] = []
    seen: set[str] = set()
    for index, label in zip(indexes, labels):
        ligand, score, source = scored[index]
        key = _ligand_stem(ligand).lower()
        if key in seen:
            continue
        seen.add(key)
        out.append((ligand, score, source, label))
    return out


def _read_consensus_ligands(
    repo_root: Path,
    run_id: str,
    pdb_id: str,
) -> list[tuple[str, float | None, str]]:
    path = run_output_dir(repo_root, "docked", run_id) / pdb_id / "consensus_docking_scores.csv"
    rows = _read_csv_rows(path)
    scored: list[tuple[str, float | None, str]] = []
    for row in rows:
        ligand = _row_ligand(row)
        score = _to_float(row.get("consensus_score"))
        if ligand and score is not None:
            scored.append((ligand, score, "consensus_score"))
    scored.sort(key=lambda item: float(item[1] or 0.0), reverse=True)
    return scored


def _read_docking_ligands(
    repo_root: Path,
    run_id: str,
    pdb_id: str,
    stage: str,
) -> list[tuple[str, float | None, str]]:
    path = run_output_dir(repo_root, "docked", run_id) / pdb_id / "docking_score_long.csv"
    rows = _read_csv_rows(path)
    scored: list[tuple[str, float | None, str]] = []
    for row in rows:
        if str(row.get("stage") or "") != str(stage):
            continue
        ligand = _row_ligand(row)
        score = _to_float(row.get("score"))
        if ligand and score is not None:
            scored.append((ligand, score, "vina_score"))
    scored.sort(key=lambda item: float(item[1] or 0.0))
    return scored


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _row_ligand(row: Mapping[str, Any]) -> str:
    for key in ("ligand", "ligand_base", "ligand_display", "drug_id"):
        value = str(row.get(key) or "").strip()
        if value:
            return value
    return ""


def _to_float(value: Any) -> float | None:
    try:
        if value is None or str(value).strip() == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _pose_path(repo_root: Path, run_id: str, pdb_id: str, ligand: str, stage: str) -> Path:
    ligand_stem = _ligand_stem(ligand)
    return (
        run_output_dir(repo_root, "docked", run_id)
        / pdb_id
        / stage
        / f"{ligand_stem}_{stage}.pdbqt"
    )


def _restore_pose_from_archive(repo_root: Path, args: argparse.Namespace, pose: Path) -> None:
    index_path = run_output_dir(repo_root, "docked", args.run_id) / "_artifact_archives" / "archive_index.json"
    if not index_path.exists():
        return
    group_id = _archive_group_for_pose(index_path, pose, args)
    if not group_id:
        return
    env = dict(os.environ)
    tmp_root = Path.home() / "atlas_tmp"
    tmp_root.mkdir(parents=True, exist_ok=True)
    env["TMPDIR"] = str(tmp_root)
    cmd = [
        sys.executable,
        "-m",
        "post_docking.artifact_retention",
        "--restore",
        "--run-id",
        str(args.run_id),
        "--restore-group",
        group_id,
        "--restore-path-glob",
        str(pose),
        "--repo-root",
        str(repo_root),
    ]
    subprocess.run(cmd, check=False, env=env)


def _archive_group_for_pose(
    index_path: Path,
    pose: Path,
    args: argparse.Namespace,
) -> str:
    import json

    payload = json.loads(index_path.read_text(encoding="utf-8"))
    wanted = str(pose.resolve())
    for group in payload.get("groups", []):
        if not _group_matches(group.get("group") or {}, args):
            continue
        entries = group.get("entries")
        if not isinstance(entries, list):
            continue
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            if str(entry.get("original_path") or "") == wanted:
                return str(group.get("group_id") or "")
    return ""


def _group_matches(group: Mapping[str, Any], args: argparse.Namespace) -> bool:
    variant = str(args.variant or "").strip().lower()
    ph = str(args.ph or "").strip().lower()
    if variant and str(group.get("variant") or "").strip().lower() != variant:
        return False
    if ph and str(group.get("ph") or "").strip().lower() != ph:
        return False
    return True


def _make_selection(
    *,
    pdb_id: str,
    ligand: str,
    score: float | None,
    score_source: str,
    pose_path: Path,
    receptor_path: Path,
    out_dir: Path,
    without_ligand: bool,
    selection_label: str,
) -> ScreenshotSelection:
    ligand_stem = _safe_name(_ligand_stem(ligand))
    suffix = "_no_ligand" if without_ligand else ""
    output_prefix = out_dir / pdb_id / ligand_stem / f"{ligand_stem}{suffix}"
    return ScreenshotSelection(
        pdb_id=pdb_id,
        ligand=ligand,
        score=score,
        score_source=score_source,
        pose_path=pose_path.resolve(),
        receptor_path=receptor_path.resolve(),
        output_prefix=output_prefix.resolve(),
        selection_label=selection_label,
    )


def _write_gallery_html(
    selections: Sequence[ScreenshotSelection],
    args: argparse.Namespace,
    out_dir: Path,
) -> Path:
    gallery_path = out_dir / "gallery.html"
    rows = "\n".join(_gallery_card(item) for item in selections)
    title = f"Atlas Screenshot Gallery: {html.escape(str(args.run_id))}"
    gallery_path.parent.mkdir(parents=True, exist_ok=True)
    gallery_path.write_text(
        f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <title>{title}</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, \"Segoe UI\", sans-serif; margin: 24px; color: #202124; }}
    h1 {{ font-size: 22px; margin-bottom: 6px; }}
    .meta {{ color: #5f6368; margin-bottom: 18px; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 18px; }}
    .card {{ border: 1px solid #dadce0; border-radius: 8px; padding: 12px; background: #fff; }}
    .card h2 {{ font-size: 16px; margin: 0 0 8px; }}
    .score {{ color: #5f6368; font-size: 13px; margin-bottom: 10px; }}
    .views {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 8px; }}
    figure {{ margin: 0; }}
    img {{ width: 100%; height: auto; border: 1px solid #eee; background: #fff; }}
    figcaption {{ font-size: 12px; color: #5f6368; margin-top: 3px; }}
    code {{ font-size: 12px; }}
  </style>
</head>
<body>
  <h1>{title}</h1>
  <div class=\"meta\">stage={html.escape(str(args.stage))} protein_style={html.escape(str(args.protein_style))} view_context={html.escape(str(args.view_context))}</div>
  <div class=\"grid\">
{rows}
  </div>
</body>
</html>
""",
        encoding="utf-8",
    )
    return gallery_path


def _gallery_card(item: ScreenshotSelection) -> str:
    score = "" if item.score is None else f"{item.score:g}"
    label = html.escape(item.selection_label or "selected")
    ligand = html.escape(item.ligand)
    pdb_id = html.escape(item.pdb_id)
    score_source = html.escape(item.score_source)
    views = []
    for view in ("front", "side", "top"):
        image = Path(f"{item.output_prefix}_{view}.png")
        data_uri = _image_data_uri(image)
        src = html.escape(data_uri)
        views.append(
            f"<figure><a href=\"{src}\"><img src=\"{src}\" alt=\"{pdb_id} {ligand} {view}\"></a>"
            f"<figcaption>{view}</figcaption></figure>"
        )
    return (
        "    <section class=\"card\">\n"
        f"      <h2>{pdb_id} / {ligand}</h2>\n"
        f"      <div class=\"score\">{label} | {score_source}"
        f"{' | score=' + html.escape(score) if score else ''}</div>\n"
        f"      <div class=\"views\">{''.join(views)}</div>\n"
        f"      <p><code>{html.escape(str(item.pose_path))}</code></p>\n"
        "    </section>"
    )


def _image_data_uri(path: Path) -> str:
    if not path.exists():
        return ""
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{payload}"


def _ligand_stem(value: str) -> str:
    name = Path(str(value or "").strip()).name
    lowered = name.lower()
    for suffix in (".pdbqt", ".pdb", ".sdf", ".mol2"):
        if lowered.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _safe_name(value: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value)
    return safe.strip("._") or "ligand"


def _parse_viewport(value: str) -> tuple[int, int]:
    raw = str(value or "").lower().strip()
    if "x" in raw:
        left, right = raw.split("x", 1)
        return max(1, int(left)), max(1, int(right))
    size = max(1, int(raw or "1600"))
    return size, size


def _print_selections(selections: Sequence[ScreenshotSelection]) -> None:
    for item in selections:
        score = "" if item.score is None else f" score={item.score:g}"
        print(
            f"{item.pdb_id}\t{item.ligand}\tpose={item.pose_path}\t"
            f"receptor={item.receptor_path}\tout={item.output_prefix}{score}"
        )
