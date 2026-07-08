from __future__ import annotations

import json
from pathlib import Path


def _notebook_source() -> str:
    notebook = Path(__file__).resolve().parents[2] / "notebooks" / "atlas_ml_exploration.ipynb"
    payload = json.loads(notebook.read_text(encoding="utf-8"))
    return "\n".join("".join(cell.get("source", [])) for cell in payload["cells"])


def test_ml_exploration_notebook_uses_run_scoped_outputs_with_legacy_fallback() -> None:
    source = _notebook_source()

    assert "REPO_ROOT" in source
    assert "outputs" in source
    assert "data" in source
    assert "RUN_ID" in source
    assert "ml" in source
    assert "training_pass" in source
    assert "LEGACY_RUN" in source
    assert "data/pilotstudy/ml_training_pass_current_defaults" in source
    assert "REPO_ROOT / 'outputs' / 'data' / RUN_ID / 'ml'" in source
