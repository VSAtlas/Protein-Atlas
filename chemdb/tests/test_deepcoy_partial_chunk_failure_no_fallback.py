from pathlib import Path
import subprocess
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
DEEPCOY_ROOT = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_ROOT) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_ROOT))

import generate_dud_library as deepcoy_lib  # type: ignore  # noqa: E402


def _write_actives(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def test_deepcoy_partial_chunk_failure_no_fallback(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("DEEPCOY_DECOYS_PER_ACTIVE", "1")
    actives_smi_path = tmp_path / "actives.smi"
    _write_actives(actives_smi_path, ["CC\tA", "CCC\tB"])
    output_dir = tmp_path / "output"
    artifact_dir = tmp_path / "artifact"
    deepcoy_run_sh = tmp_path / "deepcoy_run.sh"
    deepcoy_sdf_sh = tmp_path / "deepcoy_sdf.sh"
    config_path = tmp_path / "config.json"
    config_path.write_text("{}")
    run_log_path = tmp_path / "deepcoy_run.log"

    def _fake_run(cmd, **_kwargs):
        output_path = Path(cmd[2])
        if "chunk_000" in str(output_path):
            output_path.mkdir(parents=True, exist_ok=True)
            (output_path / "deepcoy_decoys.smi").write_text("C\tDECOY0\n")
            return subprocess.CompletedProcess(cmd, 0, stdout="ok", stderr="")
        if "chunk_001" in str(output_path):
            raise subprocess.CalledProcessError(4, cmd, output="bad", stderr="rejected")
        raise AssertionError(f"unexpected deepcoy output dir: {output_path}")

    monkeypatch.setattr(deepcoy_lib.subprocess, "run", _fake_run)

    deepcoy_lib.run_deepcoy_workflow(
        actives_smi_path,
        output_dir,
        artifact_dir,
        deepcoy_run_sh,
        deepcoy_sdf_sh,
        decoys_per_active=1,
        deepcoy_python="python",
        config_path=config_path,
        skip_sdf=True,
        ensure_sdf=False,
        restrict_data=0,
        run_log_path=run_log_path,
        deepcoy_workers=2,
        deepcoy_chunk_size=1,
        deepcoy_threads=1,
        deepcoy_keep_chunks=False,
    )

    merged_path = output_dir / "deepcoy_decoys.smi"
    assert merged_path.is_file()
    merged_text = merged_path.read_text()
    assert "DECOY0" in merged_text
    assert merged_text.strip()
    run_log_text = run_log_path.read_text()
    assert "skipping chunk" in run_log_text
