from __future__ import annotations

import json
import shutil
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

from analysis.external.pubchem_bioassay import properties_for_cids
from analysis.ml.chemical_clusters import add_butina_chemical_clusters


UPPER_BOUND_RELATIONS = {">", ">=", ">>"}
DEFAULT_FAMILIES = ("Protease", "Ion Channel", "Kinase")


def stage_target_positive_additions(
    *,
    candidates_path: str | Path,
    out_dir: str | Path,
    library_name: str,
    families: Sequence[str] = DEFAULT_FAMILIES,
    allow_non_fda: bool = False,
    resolve_pubchem: bool = True,
    retry_buffer: int = 2,
    target_cap: int = 6,
    active_nm: float = 1000.0,
    prepare: bool = True,
    force: bool = False,
    config_path: str | Path = "config.txt",
) -> dict[str, Any]:
    """Select, resolve, and prepare a sparse measured-positive add-on cohort.

    Selection is label-guided data collection, not model fitting. Candidate
    provenance and FDA/probe domain remain explicit, and the command never
    mutates the primary model-ready table.
    """

    candidates_file = Path(candidates_path)
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(candidates_file, low_memory=False)
    resolved, resolution = _resolve_structures(
        raw,
        cache_path=output / "pubchem_cid_properties.csv",
        fetch=resolve_pubchem,
    )
    eligible, cluster_summary = _eligible_candidates(
        resolved,
        families=families,
        allow_non_fda=allow_non_fda,
        active_nm=active_nm,
    )
    selected = _select_gap_fill(
        eligible,
        retry_buffer=retry_buffer,
        target_cap=target_cap,
    )
    selected = _finalize_selected(selected, library_name=library_name)

    resolved_path = output / "resolved_measured_positive_candidates.csv"
    eligible_path = output / "eligible_measured_positive_candidates.csv"
    selected_path = output / "selected_target_positive_additions.csv"
    ligand_path = output / "selected_unique_ligands.csv"
    failures_path = output / "sdf_generation_failures.csv"
    sdf_path = output / f"{library_name}.sdf"

    resolved.to_csv(resolved_path, index=False)
    eligible.to_csv(eligible_path, index=False)
    selected.to_csv(selected_path, index=False)
    _unique_ligands(selected).to_csv(ligand_path, index=False)

    failures = _write_sdf(selected, sdf_path)
    failures.to_csv(failures_path, index=False)

    prep_manifest: dict[str, Any] = {
        "prepared": False,
        "library_root": "",
        "per_pdb_libraries": 0,
        "pdbqt_available_rows": 0,
    }
    if prepare and not selected.empty and sdf_path.exists():
        prep_manifest = _prepare_libraries(
            selected,
            sdf_path=sdf_path,
            output=output,
            library_name=library_name,
            force=force,
            config_path=Path(config_path),
        )
        availability_path = Path(str(prep_manifest["availability_manifest"]))
        availability = pd.read_csv(availability_path, low_memory=False)
        selected = selected.merge(
            availability,
            on=["pdb_id", "target_gene", "ligand_base", "ligand_file_stem"],
            how="left",
            validate="one_to_one",
        )
        selected.to_csv(selected_path, index=False)


    manifest = {
        "candidates": str(candidates_file),
        "out_dir": str(output),
        "library_name": library_name,
        "families": list(families),
        "allow_non_fda": bool(allow_non_fda),
        "resolve_pubchem": bool(resolve_pubchem),
        "retry_buffer": int(retry_buffer),
        "target_cap": int(target_cap),
        "active_nm": float(active_nm),
        "candidate_rows": int(len(raw)),
        "resolved_structure_rows": int(_nonempty(resolved.get("smiles")).sum()),
        "eligible_rows": int(len(eligible)),
        "selected_pair_rows": int(len(selected)),
        "selected_unique_ligands": int(selected.get("ligand_file_stem", pd.Series(dtype=str)).nunique()),
        "selected_by_family": _count_dict(selected, "target_family"),
        "selected_by_target": _count_dict(selected, "target_gene"),
        "selected_by_source": _count_dict(selected, "source"),
        "selected_by_domain": _count_dict(selected, "ligand_domain"),
        "pubchem_resolution": resolution,
        "chemical_clustering": cluster_summary,
        "sdf_failures": int(len(failures)),
        "preparation": prep_manifest,
        "outputs": {
            "resolved_candidates": str(resolved_path),
            "eligible_candidates": str(eligible_path),
            "selected_pairs": str(selected_path),
            "selected_unique_ligands": str(ligand_path),
            "sdf": str(sdf_path),
            "sdf_failures": str(failures_path),
        },
        "policy": {
            "primary_spd_table_mutated": False,
            "candidate_absence_is_negative": False,
            "non_fda_rows_are_primary_fda_truth": False,
            "selection_uses_labels_for_data_collection_only": True,
            "required_evidence": "exact or gene-canonicalized target match; Ki/Kd/IC50 <= 1000 nM; usable structure",
            "evaluation_requirement": "retain source/target/family/chemical-cluster holdouts and report probe sensitivity separately",
        },
    }
    manifest_path = output / "target_positive_addition_install_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return manifest


def _resolve_structures(
    frame: pd.DataFrame,
    *,
    cache_path: Path,
    fetch: bool,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    work = frame.copy()
    if "pubchem_cid" not in work.columns:
        work["pubchem_cid"] = pd.Series(pd.NA, index=work.index, dtype="Int64")
    work["pubchem_cid"] = pd.to_numeric(work["pubchem_cid"], errors="coerce").astype("Int64")

    cached = pd.DataFrame(columns=["pubchem_cid", "pubchem_title", "inchikey", "smiles"])
    if cache_path.exists():
        cached = pd.read_csv(cache_path, low_memory=False)
        cached["pubchem_cid"] = pd.to_numeric(cached["pubchem_cid"], errors="coerce").astype("Int64")

    wanted = sorted(set(work.loc[work["pubchem_cid"].notna(), "pubchem_cid"].astype(int)))
    cached_ids = set(cached["pubchem_cid"].dropna().astype(int))
    cached_complete_ids = set(
        cached.loc[_nonempty(cached.get("smiles")), "pubchem_cid"].dropna().astype(int)
    )
    missing = sorted(set(wanted) - cached_complete_ids)
    fetched = properties_for_cids(missing) if fetch and missing else pd.DataFrame()
    properties = pd.concat([cached, fetched], ignore_index=True, sort=False)
    if not properties.empty:
        properties["pubchem_cid"] = pd.to_numeric(properties["pubchem_cid"], errors="coerce").astype("Int64")
        properties = properties.dropna(subset=["pubchem_cid"]).drop_duplicates("pubchem_cid", keep="last")
        properties.to_csv(cache_path, index=False)

    if properties.empty:
        work["structure_resolution_status"] = _nonempty(work.get("smiles")).map(
            {True: "source_structure", False: "missing_structure"}
        )
        return work, {
            "requested_cids": len(wanted),
            "cached_cids": len(cached_ids),
            "fetched_cids": 0,
            "resolved_cids": 0,
        }

    props = properties.rename(
        columns={
            "pubchem_title": "_pubchem_title",
            "inchikey": "_pubchem_inchikey",
            "smiles": "_pubchem_smiles",
        }
    )
    work = work.merge(props, on="pubchem_cid", how="left")

    source_key = work.get("inchikey", pd.Series("", index=work.index)).map(_clean)
    pubchem_key = work.get("_pubchem_inchikey", pd.Series("", index=work.index)).map(_clean)
    conflict = (
        source_key.ne("")
        & pubchem_key.ne("")
        & source_key.map(_inchikey_connectivity).ne(pubchem_key.map(_inchikey_connectivity))
    )
    source_smiles = work.get("smiles", pd.Series("", index=work.index)).map(_clean)
    pubchem_smiles = work.get("_pubchem_smiles", pd.Series("", index=work.index)).map(_clean)
    work["smiles"] = source_smiles.where(source_smiles.ne(""), pubchem_smiles)
    work["inchikey"] = source_key.where(source_key.ne(""), pubchem_key)
    work["structure_resolution_status"] = "missing_structure"
    work.loc[source_smiles.ne(""), "structure_resolution_status"] = "source_structure"
    work.loc[source_smiles.eq("") & pubchem_smiles.ne(""), "structure_resolution_status"] = "pubchem_resolved"
    work.loc[conflict, "structure_resolution_status"] = "source_pubchem_conflict"

    title = work.get("_pubchem_title", pd.Series("", index=work.index)).map(_clean)
    for field in ("drug_name", "display_name", "generic_name"):
        if field not in work.columns:
            work[field] = ""
        current = work[field].map(_clean)
        work[field] = current.where(current.ne(""), title)
    work = work.drop(
        columns=["_pubchem_title", "_pubchem_inchikey", "_pubchem_smiles"],
        errors="ignore",
    )
    return work, {
        "requested_cids": len(wanted),
        "cached_cids": len(cached_ids),
        "fetched_cids": int(len(fetched)),
        "resolved_cids": int(properties["pubchem_cid"].nunique()),
        "structure_conflicts": int(conflict.sum()),
    }


def _eligible_candidates(
    frame: pd.DataFrame,
    *,
    families: Sequence[str],
    allow_non_fda: bool,
    active_nm: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    work = frame.copy()
    if families:
        work = work[work.get("target_family").isin(families)].copy()

    relation = work.get("activity_relation", pd.Series("", index=work.index)).astype(str).str.strip()
    endpoint = work.get("activity_type", pd.Series("", index=work.index)).astype(str).str.upper().str.strip()
    activity = pd.to_numeric(work.get("activity_nM"), errors="coerce")
    training_allowed = work.get("training_allowed", pd.Series(True, index=work.index)).map(_truthy)
    benchmark_only = work.get("benchmark_only", pd.Series(False, index=work.index)).map(_truthy)
    valid_structure = _nonempty(work.get("smiles"))
    no_conflict = work.get(
        "structure_resolution_status",
        pd.Series("source_structure", index=work.index),
    ).ne("source_pubchem_conflict")
    exact_target = _nonempty(work.get("target_gene")) & _nonempty(work.get("pdb_id"))
    eligible = (
        training_allowed
        & ~benchmark_only
        & activity.le(active_nm)
        & ~relation.isin(UPPER_BOUND_RELATIONS)
        & endpoint.isin({"KI", "KD", "IC50"})
        & valid_structure
        & no_conflict
        & exact_target
    )
    if not allow_non_fda:
        eligible &= work.get("production_truth_allowed", pd.Series(False, index=work.index)).map(_truthy)
    work = work[eligible].copy()
    work["_compound_key"] = work.apply(_compound_key, axis=1)
    work = work.drop_duplicates(["target_uniprot", "pdb_id", "_compound_key"], keep="first")
    work, cluster_summary = add_butina_chemical_clusters(work, distance_threshold=0.35)
    return work.drop(columns=["_compound_key"], errors="ignore").reset_index(drop=True), cluster_summary


def _select_gap_fill(
    eligible: pd.DataFrame,
    *,
    retry_buffer: int,
    target_cap: int,
) -> pd.DataFrame:
    selected_rows: list[pd.Series] = []
    if eligible.empty:
        return eligible.copy()

    for family, family_rows in eligible.groupby("target_family", dropna=False):
        gap = pd.to_numeric(family_rows.get("family_positive_gap"), errors="coerce").max()
        gap_int = max(0, int(gap)) if pd.notna(gap) else 0
        if gap_int <= 0:
            continue
        need = gap_int + max(0, int(retry_buffer))
        pool = family_rows.copy()
        target_counts: Counter[str] = Counter()
        cluster_counts: Counter[str] = Counter()

        while len(pool) and sum(target_counts.values()) < need:
            scored: list[tuple[tuple[float, ...], Any]] = []
            for index, row in pool.iterrows():
                target = _clean(row.get("target_gene"))
                cluster = _clean(row.get("chemical_cluster")) or "unassigned"
                if target_counts[target] >= max(1, int(target_cap)):
                    continue
                target_gap = _as_number(row.get("target_positive_gap"))
                target_need_bonus = 100.0 if target_counts[target] < target_gap else 0.0
                new_target_bonus = 75.0 if target_counts[target] == 0 else 0.0
                new_cluster_bonus = 20.0 if cluster_counts[cluster] == 0 else 0.0
                source_bonus = 5.0 if _clean(row.get("source")) not in {
                    _clean(item.get("source")) for item in selected_rows if _clean(item.get("target_family")) == _clean(family)
                } else 0.0
                potency = -_as_number(row.get("activity_nM"), default=active_default())
                score = (
                    target_need_bonus + new_target_bonus + new_cluster_bonus + source_bonus,
                    -float(target_counts[target]),
                    -float(cluster_counts[cluster]),
                    potency,
                    -float(index),
                )
                scored.append((score, index))
            if not scored:
                break
            _, chosen_index = max(scored, key=lambda item: item[0])
            chosen = pool.loc[chosen_index].copy()
            selected_rows.append(chosen)
            target_counts[_clean(chosen.get("target_gene"))] += 1
            cluster_counts[_clean(chosen.get("chemical_cluster")) or "unassigned"] += 1
            pool = pool.drop(index=chosen_index)

    if not selected_rows:
        return eligible.iloc[0:0].copy()
    return pd.DataFrame(selected_rows).reset_index(drop=True)


def active_default() -> float:
    return 1_000_000_000.0


def _finalize_selected(frame: pd.DataFrame, *, library_name: str) -> pd.DataFrame:
    selected = frame.copy()
    if selected.empty:
        return selected
    from prep_ligands.prep_ligands_runtime import sanitize_ligand_name_for_filename

    if "ligand_base" not in selected.columns:
        selected["ligand_base"] = selected["drug_id"]
    selected["ligand_file_stem"] = selected["ligand_base"].map(
        lambda value: sanitize_ligand_name_for_filename(_clean(value))
    )
    selected["library_name"] = library_name
    selected["primary_fda_claim_allowed"] = selected.get(
        "production_truth_allowed",
        pd.Series(False, index=selected.index),
    ).map(_truthy)
    selected["probe_sensitivity_only"] = ~selected["primary_fda_claim_allowed"]
    selected["selection_reason"] = (
        selected.get("selection_reason", pd.Series("", index=selected.index)).map(_clean)
        + ";automated_family_gap_fill;chemical_diversity_aware"
    ).str.strip(";")
    selected["training_domain"] = selected["training_domain"].where(
        selected["primary_fda_claim_allowed"],
        "external_probe_sensitivity",
    )
    return selected


def _unique_ligands(selected: pd.DataFrame) -> pd.DataFrame:
    fields = [
        "ligand_base",
        "ligand_file_stem",
        "display_name",
        "generic_name",
        "pubchem_cid",
        "inchikey",
        "smiles",
        "ligand_domain",
        "training_domain",
        "primary_fda_claim_allowed",
        "probe_sensitivity_only",
    ]
    columns = [field for field in fields if field in selected.columns]
    if not columns:
        return pd.DataFrame()
    return selected[columns].drop_duplicates("ligand_file_stem").reset_index(drop=True)


def _write_sdf(selected: pd.DataFrame, sdf_path: Path) -> pd.DataFrame:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    failures: list[dict[str, str]] = []
    unique = _unique_ligands(selected)
    writer = Chem.SDWriter(str(sdf_path))
    try:
        for _, row in unique.iterrows():
            stem = _clean(row.get("ligand_file_stem"))
            smiles = _clean(row.get("smiles"))
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                failures.append({"ligand_file_stem": stem, "reason": "invalid_smiles"})
                continue
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            status = AllChem.EmbedMolecule(mol, params)
            if status != 0:
                params.useRandomCoords = True
                status = AllChem.EmbedMolecule(mol, params)
            if status != 0:
                failures.append({"ligand_file_stem": stem, "reason": "rdkit_embedding_failed"})
                continue
            try:
                if AllChem.MMFFHasAllMoleculeParams(mol):
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
                else:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                pass
            mol.SetProp("_Name", stem)
            for field in ("ligand_base", "display_name", "pubchem_cid", "inchikey", "ligand_domain"):
                value = _clean(row.get(field))
                if value:
                    mol.SetProp(f"ATLAS_{field.upper()}", value)
            writer.write(mol)
    finally:
        writer.close()
    return pd.DataFrame(failures, columns=["ligand_file_stem", "reason"])


def _prepare_libraries(
    selected: pd.DataFrame,
    *,
    sdf_path: Path,
    output: Path,
    library_name: str,
    force: bool,
    config_path: Path,
) -> dict[str, Any]:
    from config.runtime_config import load_config
    from path_router import make_paths
    from prep_ligands.library_index import LibraryIndex
    from prep_ligands.prep_ligands_bulk_meeko import prep_ligands_with_meeko

    repo_root = Path(__file__).resolve().parents[2]
    config = config_path if config_path.is_absolute() else repo_root / config_path
    cfg = load_config(str(config), base_dir=repo_root)
    first_pdb = _clean(selected.iloc[0]["pdb_id"]).upper()
    paths = make_paths(cfg, first_pdb, f"{first_pdb}.pdb")
    all_library = paths.prepped_root / f"{library_name}_all"
    all_library.mkdir(parents=True, exist_ok=True)

    prep_ligands_with_meeko(
        force=force,
        in_sdf=sdf_path,
        mol2_dir=output / "ligand_intermediates",
        out_pdbqt_dir=all_library,
        status_log=output / "ligand_prep_status.tsv",
    )
    all_index = LibraryIndex()
    relative_pdbqts = sorted(
        path.relative_to(all_library)
        for path in all_library.rglob("*.pdbqt")
        if path.is_file()
    )
    all_index.write_manifest_for_root(all_library, relative_pdbqts)
    all_index.load([all_library])

    copy_rows: list[dict[str, Any]] = []
    availability_rows: list[dict[str, Any]] = []
    target_dirs: set[Path] = set()
    library_map: dict[str, str] = {}
    for _, row in selected.iterrows():
        pdb_id = _clean(row.get("pdb_id")).upper()
        stem = _clean(row.get("ligand_file_stem"))
        source = all_index.lookup(stem, [all_library])
        target_dir = paths.prepped_root / f"{library_name}_{pdb_id}"
        target_dir.mkdir(parents=True, exist_ok=True)
        target_dirs.add(target_dir)
        library_map[pdb_id] = target_dir.name
        destination = target_dir / f"{stem}.pdbqt"
        available = source is not None and source.is_file()
        if available:
            shutil.copy2(source, destination)
            sidecar = source.with_suffix(".ligprep_source.json")
            if sidecar.is_file():
                shutil.copy2(sidecar, destination.with_suffix(".ligprep_source.json"))
        copy_rows.append(
            {
                "pdb_id": pdb_id,
                "library": target_dir.name,
                "ligand_base": row.get("ligand_base"),
                "ligand_file_stem": stem,
                "source": str(source) if source else "",
                "destination": str(destination),
                "copied": bool(available),
            }
        )
        availability_rows.append(
            {
                "pdb_id": pdb_id,
                "target_gene": row.get("target_gene"),
                "ligand_base": row.get("ligand_base"),
                "ligand_file_stem": stem,
                "pdbqt_path": str(destination),
                "pdbqt_available": bool(available and destination.is_file()),
            }
        )

    for target_dir in target_dirs:
        index = LibraryIndex()
        index.load([target_dir])

    copies = pd.DataFrame(copy_rows)
    availability = pd.DataFrame(availability_rows)
    copies.to_csv(output / "per_pdb_library_copy_manifest.csv", index=False)
    availability.to_csv(output / "pdbqt_availability_after_prep.csv", index=False)
    library_map_path = output / "per_pdb_library_map.json"
    library_map_path.write_text(json.dumps(library_map, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "prepared": True,
        "library_root": str(all_library),
        "per_pdb_libraries": len(target_dirs),
        "pdbqt_available_rows": int(availability.get("pdbqt_available", pd.Series(dtype=bool)).sum()),
        "pdbqt_missing_rows": int((~availability.get("pdbqt_available", pd.Series(dtype=bool))).sum()),
        "copy_manifest": str(output / "per_pdb_library_copy_manifest.csv"),
        "availability_manifest": str(output / "pdbqt_availability_after_prep.csv"),
        "library_map": str(library_map_path),
    }


def _compound_key(row: pd.Series) -> str:
    inchikey = _inchikey_connectivity(_clean(row.get("inchikey")))
    if inchikey:
        return f"inchikey:{inchikey}"
    return f"smiles:{_clean(row.get('smiles'))}"


def _inchikey_connectivity(value: str) -> str:
    return value.split("-", 1)[0].upper() if value else ""


def _clean(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _nonempty(series: pd.Series | None) -> pd.Series:
    if series is None:
        return pd.Series(dtype=bool)
    return series.map(_clean).ne("")


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return _clean(value).lower() in {"1", "true", "yes", "y"}


def _as_number(value: Any, *, default: float = 0.0) -> float:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return float(numeric) if pd.notna(numeric) else float(default)


def _count_dict(frame: pd.DataFrame, column: str) -> dict[str, int]:
    if frame.empty or column not in frame.columns:
        return {}
    return {
        str(key): int(value)
        for key, value in frame.groupby(column, dropna=False).size().to_dict().items()
    }
