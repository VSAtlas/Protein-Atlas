import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chemdb.target_difficulty import bucket_from_auc, difficulty_from_metrics


def test_bucket_from_auc_thresholds():
    assert bucket_from_auc(0.569) == "hard"
    assert bucket_from_auc(0.641) == "medium"
    assert bucket_from_auc(0.819) == "easy"
    assert bucket_from_auc(float("nan")) == "degenerate"


def test_difficulty_from_metrics_series():
    metrics = pd.Series({"ROC_AUC": 0.641, "N": 123, "n_actives": 45})
    td = difficulty_from_metrics("1BCD", metrics)
    assert td.difficulty == "medium"
    assert math.isclose(td.roc_auc, 0.641, rel_tol=1e-6)
    assert td.N == 123
    assert td.n_actives == 45


def test_difficulty_from_missing_metrics():
    td = difficulty_from_metrics("X000", None)
    assert td.difficulty == "degenerate"
    assert math.isnan(td.roc_auc)
    assert td.N == 0
    assert td.n_actives == 0
