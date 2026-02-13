import sys
from pathlib import Path

DEEPCOY_DIR = Path(__file__).resolve().parents[2] / "DeepCoy_duds"
if str(DEEPCOY_DIR) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_DIR))

from DeepCoy_duds.generate_dud_library import merge_deepcoy_outputs, prepare_deepcoy_chunks  # noqa: E402


def _write_actives(tmp_path: Path, lines):
    actives_path = tmp_path / "actives.smi"
    actives_path.write_text("".join(lines))
    return actives_path


def test_prepare_deepcoy_chunks_creates_expected_chunks(tmp_path):
    lines = ["SMI1\tA\n", "SMI2\tB\n", "SMI3\tC\n", "SMI4\tD\n", "SMI5\tE\n"]
    actives_path = _write_actives(tmp_path, lines)
    output_dir = tmp_path / "out" / "target"
    artifact_dir = tmp_path / "artifacts"
    chunk_output_root = output_dir.parent / ".deepcoy_chunks" / output_dir.name
    chunk_artifact_root = artifact_dir / "chunks"

    actives_lines, chunk_plans = prepare_deepcoy_chunks(
        actives_path,
        output_dir,
        artifact_dir,
        chunk_output_root,
        chunk_artifact_root,
        chunk_size=2,
        use_chunking=True,
    )

    assert actives_lines == lines
    assert len(chunk_plans) == 3
    assert chunk_plans[0].actives_path.read_text().splitlines() == [
        line.rstrip("\n") for line in lines[:2]
    ]
    assert chunk_plans[1].actives_path.read_text().splitlines() == [
        line.rstrip("\n") for line in lines[2:4]
    ]
    assert chunk_plans[2].actives_path.read_text().splitlines() == [
        line.rstrip("\n") for line in lines[4:]
    ]


def test_merge_deepcoy_outputs_deduplicates_and_orders(tmp_path):
    actives_path = _write_actives(tmp_path, ["A\n", "B\n"])
    output_dir = tmp_path / "out" / "target"
    artifact_dir = tmp_path / "artifacts"
    chunk_output_root = output_dir.parent / ".deepcoy_chunks" / output_dir.name
    chunk_artifact_root = artifact_dir / "chunks"
    _, chunk_plans = prepare_deepcoy_chunks(
        actives_path,
        output_dir,
        artifact_dir,
        chunk_output_root,
        chunk_artifact_root,
        chunk_size=1,
        use_chunking=True,
    )

    decoy_chunks = [
        ["D1\tchunk0\n", "D2\tchunk0\n"],
        ["D2\tchunk1\n", "D3\tchunk1\n"],
    ]
    for plan, decoy_lines in zip(chunk_plans, decoy_chunks):
        decoy_path = plan.output_dir / "deepcoy_decoys.smi"
        plan.output_dir.mkdir(parents=True, exist_ok=True)
        decoy_path.write_text("".join(decoy_lines))

    merged_path = output_dir / "deepcoy_decoys.smi"
    produced = merge_deepcoy_outputs(
        chunk_plans,
        merged_path,
        decoys_per_active=2,
        total_actives=2,
        keep_chunks=False,
        chunk_output_root=chunk_output_root,
        chunk_artifact_root=chunk_artifact_root,
    )

    assert produced == 3
    assert merged_path.read_text().splitlines() == [
        "D1\tchunk0",
        "D2\tchunk0",
        "D3\tchunk1",
    ]
    assert not chunk_output_root.exists()
    assert not chunk_artifact_root.exists()


def test_prepare_chunks_workers_one_matches_single_run(tmp_path):
    lines = ["AA\tfirst\n", "BB\tsecond\n"]
    actives_path = _write_actives(tmp_path, lines)
    output_dir = tmp_path / "out" / "target"
    artifact_dir = tmp_path / "artifacts"
    chunk_output_root = output_dir.parent / ".deepcoy_chunks" / output_dir.name
    chunk_artifact_root = artifact_dir / "chunks"

    actives_lines, plans = prepare_deepcoy_chunks(
        actives_path,
        output_dir,
        artifact_dir,
        chunk_output_root,
        chunk_artifact_root,
        chunk_size=3,
        use_chunking=False,
    )

    assert actives_lines == lines
    assert len(plans) == 1
    plan = plans[0]
    assert plan.actives_path == actives_path
    assert plan.output_dir == output_dir
    assert plan.artifact_dir == artifact_dir
    assert plan.line_count == len(lines)
    assert plan.is_chunk is False
    assert not chunk_output_root.exists()
    assert not chunk_artifact_root.exists()
