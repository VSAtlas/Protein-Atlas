from pathlib import Path

from analysis.reporting.ligand_annotation_groups import (
    build_ligand_grouping_meta,
    resolve_ligand_annotations,
)


def test_resolve_ligand_annotations_override_hits(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    (repo_root / "chemdb").mkdir(parents=True)

    levonorgestrel = resolve_ligand_annotations(
        repo_root, ligand_name="levonorgestrel", ligand_base="levonorgestrel"
    )
    assert levonorgestrel["chemotype_primary"] == "Steroid-like"
    assert levonorgestrel["drug_effect_primary"] == "Steroid hormone / endocrine modulator"
    assert "steroidal core" in levonorgestrel["motif_tags"]
    assert levonorgestrel["ligand_annotation_confidence"] == "high"

    itraconazole = resolve_ligand_annotations(
        repo_root, ligand_name="itraconazole", ligand_base="itraconazole"
    )
    assert itraconazole["chemotype_primary"] == "Azole-like"
    assert itraconazole["drug_effect_primary"] == "Antifungal"
    assert "azole ring" in itraconazole["motif_tags"]


def test_build_ligand_grouping_meta_counts_blocks(tmp_path: Path) -> None:
    repo_root = tmp_path / "repo"
    (repo_root / "chemdb").mkdir(parents=True)
    rows = [
        {
            "ligand_name": "triamcinolone",
            "ligand_display": "triamcinolone",
            "ligand_base": "triamcinolone",
            "library": "FDA",
        },
        {
            "ligand_name": "itraconazole",
            "ligand_display": "itraconazole",
            "ligand_base": "itraconazole",
            "library": "FDA",
        },
    ]
    meta, summary = build_ligand_grouping_meta(
        repo_root, rows, ["triamcinolone", "itraconazole"]
    )
    assert meta["triamcinolone"]["drug_effect_primary"] == "Glucocorticoid"
    assert meta["itraconazole"]["chemotype_primary"] == "Azole-like"
    chemotype_labels = [item["label"] for item in summary["chemotype_blocks"]]
    effect_labels = [item["label"] for item in summary["drug_effect_blocks"]]
    assert "Corticosteroid-like" in chemotype_labels
    assert "Azole-like" in chemotype_labels
    assert "Glucocorticoid" in effect_labels
    assert "Antifungal" in effect_labels
