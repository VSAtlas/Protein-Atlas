from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
DEEPCOY_ROOT = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_ROOT) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_ROOT))

from DeepCoy_duds.generate_dud_library import (  # type: ignore  # noqa: E402
    prepare_deepcoy_chunks,
    resolve_deepcoy_chunk_roots,
)


def _write_actives(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def test_resolve_chunk_roots_keep_chunks_preserves_existing_runs(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setenv("DEEPCOY_DECOYS_PER_ACTIVE", "1")
    out_root = tmp_path / "out"
    artifact_dir = tmp_path / "artifact"
    label = "label"
    actives_path = out_root / label / f"{label}_actives.smi"
    _write_actives(actives_path, ["CC\tA", "CCC\tB"])

    old_chunk_marker = out_root / ".deepcoy_chunks" / label / "oldrun" / "marker.txt"
    old_artifact_marker = artifact_dir / "chunks" / "oldrun" / "marker.txt"
    old_chunk_marker.parent.mkdir(parents=True, exist_ok=True)
    old_chunk_marker.write_text("keep")
    old_artifact_marker.parent.mkdir(parents=True, exist_ok=True)
    old_artifact_marker.write_text("keep")

    run_tag = "testrun"
    (
        chunk_output_root,
        chunk_artifact_root,
        resolved_label,
    ) = resolve_deepcoy_chunk_roots(
        out_root / label,
        artifact_dir,
        actives_path,
        keep_chunks=True,
        run_tag=run_tag,
    )

    assert resolved_label == label
    assert chunk_output_root == out_root / ".deepcoy_chunks" / label / run_tag
    assert chunk_artifact_root == artifact_dir / "chunks" / run_tag

    _, chunk_plans = prepare_deepcoy_chunks(
        actives_path,
        out_root / label,
        artifact_dir,
        chunk_output_root,
        chunk_artifact_root,
        chunk_size=1,
        use_chunking=True,
    )
    assert chunk_plans, "chunk plans should be produced for chunking path"
    assert old_chunk_marker.exists()
    assert old_artifact_marker.exists()


def test_resolve_chunk_roots_legacy_paths_when_not_keeping_chunks(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setenv("DEEPCOY_DECOYS_PER_ACTIVE", "1")
    out_root = tmp_path / "out"
    artifact_dir = tmp_path / "artifact"
    label = "label"
    actives_path = out_root / label / f"{label}_actives.smi"
    _write_actives(actives_path, ["CC\tA"])

    (
        chunk_output_root,
        chunk_artifact_root,
        resolved_label,
    ) = resolve_deepcoy_chunk_roots(
        out_root / label,
        artifact_dir,
        actives_path,
        keep_chunks=False,
        run_tag="ignored",
    )

    assert resolved_label == label
    assert chunk_output_root == out_root / ".deepcoy_chunks" / label
    assert chunk_artifact_root == artifact_dir / "chunks"
