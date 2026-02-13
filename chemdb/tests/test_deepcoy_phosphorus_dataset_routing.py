from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
DEEPCOY_ROOT = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_ROOT) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_ROOT))

from DeepCoy_duds import generate_dud_library as deepcoy_lib  # type: ignore  # noqa: E402


def _write_actives(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def test_deepcoy_phosphorus_dataset_routing(tmp_path: Path) -> None:
    actives_path = tmp_path / "actives.smi"
    _write_actives(actives_path, ["CP\tP1", "CC\tC1"])
    output_dir = tmp_path / "out"
    artifact_dir = tmp_path / "artifact"

    chunk_output_root, chunk_artifact_root, _ = deepcoy_lib.resolve_deepcoy_chunk_roots(
        output_dir,
        artifact_dir,
        actives_path,
        keep_chunks=False,
        run_tag=None,
    )
    _, chunk_plans = deepcoy_lib.prepare_deepcoy_chunks(
        actives_path,
        output_dir,
        artifact_dir,
        chunk_output_root,
        chunk_artifact_root,
        chunk_size=1,
        use_chunking=True,
    )

    default_model = tmp_path / "default_model.pickle"
    phosphorus_model = tmp_path / "phosphorus_model.pickle"
    default_model.write_text("x")
    phosphorus_model.write_text("x")

    deepcoy_lib.assign_chunk_metadata(
        chunk_plans,
        phosphorus_model,
        default_model,
        base_seed=None,
        seed_per_chunk=False,
        run_log_path=None,
    )

    plan_phosphorus = next(plan for plan in chunk_plans if plan.index == 0)
    plan_default = next(plan for plan in chunk_plans if plan.index == 1)

    assert plan_phosphorus.dataset == "zinc_phosphorus"
    assert plan_phosphorus.model_path == phosphorus_model
    assert plan_default.dataset == "zinc"
    assert plan_default.model_path == default_model
