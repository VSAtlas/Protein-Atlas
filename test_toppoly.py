# smoke_test.py
import math
from typing import Dict, Tuple, Any, List

def score_key(item: Tuple[str, Dict[str, Any]]) -> float:
    _, rec = item
    s = rec.get("score")
    if s is None:
        return math.inf
    return float(s) if rec.get("valid", False) else float(s) + 1e-6

def select_top_polypharma(final_scores: Dict[str, Dict[str, Any]], k: int = 20) -> List[Tuple[str, Dict[str, Any]]]:
    return sorted(final_scores.items(), key=score_key)[:k]

def main():
    score_history = {
        "stage5": {
            "lig1": {"score": -8.3, "valid": True},
            "lig2": {"score": -6.0, "valid": False},
            "lig3": {"score": None, "valid": False},
            "lig4": {"score": -8.3, "valid": False},  # tie vs lig1, invalid should come after
        }
    }
    final_stage = "stage5"
    top = select_top_polypharma(score_history[final_stage], k=3)
    print("TOP:", top)
    assert top[0][0] == "lig1"
    assert top[1][0] in {"lig2", "lig4"}
    assert top[-1][0] != "lig3" or top[-1][1]["score"] is None

if __name__ == "__main__":
    main()
