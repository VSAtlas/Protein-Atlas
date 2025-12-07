import tempfile
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dud_eval import compute_decoy_stats_from_long_csv, guess_ligfile_col, guess_score_col


def main() -> None:
    tmpdir = Path(tempfile.mkdtemp(prefix="tscore_demo_"))
    dud_csv = tmpdir / "dud_docking_score_long.csv"
    fda_csv = tmpdir / "docking_score_long.csv"

    # Synthetic DUD: 3 decoy ligands, each with two poses
    dud_df = pd.DataFrame(
        {
            "ligand_file": [
                "decoy1_decoy_1.pdbqt",
                "decoy1_decoy_2.pdbqt",
                "decoy2_decoy_1.pdbqt",
                "decoy2_decoy_2.pdbqt",
                "decoy3_decoy_1.pdbqt",
                "decoy3_decoy_2.pdbqt",
            ],
            "score": [-7.0, -6.5, -8.0, -7.8, -9.0, -8.9],
        }
    )
    dud_df.to_csv(dud_csv, index=False)

    # Synthetic FDA: 2 ligands, each with two poses
    fda_df = pd.DataFrame(
        {
            "ligand_file": [
                "ligA_active_1.pdbqt",
                "ligA_active_2.pdbqt",
                "ligB_active_1.pdbqt",
                "ligB_active_2.pdbqt",
            ],
            "score": [-9.5, -9.2, -7.5, -7.0],
        }
    )
    fda_df.to_csv(fda_csv, index=False)

    # Compute decoy stats from the DUD file
    mu, sigma, n_decoys = compute_decoy_stats_from_long_csv(dud_csv)
    print(f"Decoy stats: n={n_decoys}, mean={mu:.3f}, std={sigma:.3f}")

    # Annotate the FDA file in the same way as annotate_fda_long_csv_with_t_scores_vs_decoys,
    # but without needing cfg/make_paths.
    df = pd.read_csv(fda_csv)
    lig_col = guess_ligfile_col(df, None)
    score_col = guess_score_col(df, None)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    best = (
        df.groupby(lig_col, as_index=False)
          .agg(best_score=(score_col, "min"))
    )
    best["t_vs_decoys"] = (mu - best["best_score"]) / sigma
    t_map = dict(zip(best[lig_col], best["t_vs_decoys"]))
    df["t_vs_decoys"] = df[lig_col].map(t_map)

    print("Annotated FDA scores:")
    print(df)


if __name__ == "__main__":
    main()
