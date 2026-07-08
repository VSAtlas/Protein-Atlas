from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.adrecs_target import load_adrecs_target
from analysis.external.ctd import load_ctd
from analysis.external.opentargets_safety import load_opentargets_safety
from analysis.external.reactome import load_reactome
from analysis.external.sider import load_sider


def _node_id(node_type: str, value: Any) -> str:
    return f"{node_type}:{str(value).strip()}"


def _json_payload(metadata: dict[str, Any] | None) -> str:
    clean = {key: value for key, value in (metadata or {}).items() if not pd.isna(value)}
    return json.dumps(clean, sort_keys=True)


def _add_node(nodes: dict[str, dict[str, object]], node_type: str, label: Any, source: str, metadata: dict[str, Any] | None = None) -> str:
    node_id = _node_id(node_type, label)
    nodes.setdefault(
        node_id,
        {"node_id": node_id, "node_type": node_type, "label": label, "source": source, "metadata_json": _json_payload(metadata)},
    )
    return node_id


def _ctd_id(value: Any) -> str | None:
    return str(value).strip() if pd.notna(value) and str(value).strip() else None


def _edge(src: str, dst: str, edge_type: str, evidence_source: str, confidence: float = 1.0, pubmed_ids: Any = "", metadata: dict[str, Any] | None = None) -> dict[str, object]:
    return {
        "source_node_id": src,
        "target_node_id": dst,
        "edge_type": edge_type,
        "evidence_source": evidence_source,
        "confidence": confidence,
        "pubmed_ids": "" if pd.isna(pubmed_ids) else pubmed_ids,
        "metadata_json": _json_payload(metadata),
    }


def build_mechanism_graph(
    pair_table_path: str | Path,
    sider_path: str | Path | None,
    ctd_path: str | Path | None,
    opentargets_safety_path: str | Path | None,
    reactome_path: str | Path | None,
    spd_path: str | Path | None,
    out_edges_path: str | Path,
    out_nodes_path: str | Path,
    adrecs_target_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pair = pd.read_csv(pair_table_path)
    nodes: dict[str, dict[str, object]] = {}
    edges: list[dict[str, object]] = []
    for _idx, row in pair.iterrows():
        drug = _add_node(nodes, "drug", row["drug_id"], "atlas_pair_table")
        target = _add_node(nodes, "target", row["target_id"], "atlas_pair_table")
        edges.append(_edge(drug, target, "drug_predicted_target", "Atlas", float(row.get("atlas_score", 1) or 1), metadata={"pdb_id": row.get("pdb_id", "")}))
        if pd.notna(row.get("pdb_id")) and str(row.get("pdb_id")):
            pdb = _add_node(nodes, "pdb", row["pdb_id"], "atlas_pair_table")
            edges.append(_edge(target, pdb, "target_structure", "Atlas"))
    if sider_path and Path(sider_path).exists():
        for _idx, row in load_sider(sider_path).iterrows():
            drug = _add_node(nodes, "drug", row["drug_id"], "SIDER")
            adr = _add_node(nodes, "adr", row["adr"], "SIDER")
            edges.append(_edge(drug, adr, "drug_has_adr", row.get("source", "SIDER")))
    if opentargets_safety_path and Path(opentargets_safety_path).exists():
        for _idx, row in load_opentargets_safety(opentargets_safety_path).iterrows():
            target = _add_node(nodes, "target", row["target_id"], "OpenTargetsSafety")
            adr = _add_node(nodes, "adr", row["adr"], "OpenTargetsSafety")
            conf = pd.to_numeric(pd.Series([row.get("confidence", 1.0)]), errors="coerce").fillna(1.0).iloc[0]
            edges.append(_edge(target, adr, "target_has_safety_evidence", row.get("source", "OpenTargetsSafety"), float(conf), row.get("pubmed_ids", "")))
    if adrecs_target_path and Path(adrecs_target_path).exists():
        for _idx, row in load_adrecs_target(adrecs_target_path).iterrows():
            target = _add_node(nodes, "target", row["target_id"], "ADReCS-Target")
            adr = _add_node(nodes, "adr", row["adr"], "ADReCS-Target")
            edges.append(
                _edge(
                    target,
                    adr,
                    "target_has_safety_evidence",
                    row.get("source", "ADReCS-Target"),
                    float(row.get("confidence", 1.0) or 1.0),
                    row.get("pubmed_ids", ""),
                    metadata={
                        "_source_badd_tid": row.get("_source_badd_tid", ""),
                        "_source_adr_id": row.get("_source_adr_id", ""),
                        "_source_adrecs_id": row.get("_source_adrecs_id", ""),
                        "_source_drug_name": row.get("_source_drug_name", ""),
                    },
                )
            )
    if reactome_path and Path(reactome_path).exists():
        for _idx, row in load_reactome(reactome_path).iterrows():
            target = _add_node(nodes, "target", row["target_id"], "Reactome")
            pathway = _add_node(nodes, "pathway", row["pathway"], "Reactome", {"pathway_id": row.get("pathway_id", "")})
            edges.append(_edge(target, pathway, "target_in_pathway", row.get("source", "Reactome")))
            if pd.notna(row.get("adr")) and str(row.get("adr")):
                adr = _add_node(nodes, "adr", row["adr"], "Reactome")
                edges.append(_edge(pathway, adr, "pathway_linked_to_adr", row.get("source", "Reactome")))
    if ctd_path and Path(ctd_path).exists():
        for _idx, row in load_ctd(ctd_path).iterrows():
            chemical_id = _ctd_id(row.get("chemical_id"))
            gene_id = _ctd_id(row.get("gene_id"))
            disease_id = _ctd_id(row.get("disease_id"))
            chemical = _add_node(nodes, "drug", chemical_id, "CTD") if chemical_id else None
            gene = _add_node(nodes, "target", gene_id, "CTD") if gene_id else None
            disease = _add_node(nodes, "disease", disease_id, "CTD") if disease_id else None
            if chemical and gene:
                edges.append(_edge(chemical, gene, "chemical_gene_interaction", row.get("source", "CTD"), pubmed_ids=row.get("pubmed_ids", "")))
            if chemical and disease:
                edges.append(_edge(chemical, disease, "chemical_disease_association", row.get("source", "CTD"), pubmed_ids=row.get("pubmed_ids", "")))
            if gene and disease:
                edges.append(_edge(gene, disease, "gene_disease_association", row.get("source", "CTD"), pubmed_ids=row.get("pubmed_ids", "")))
    if spd_path and Path(spd_path).exists():
        spd = pd.read_csv(spd_path)
        for _idx, row in spd.iterrows():
            drug = _add_node(nodes, "drug", row["drug_id"], "SPD")
            target = _add_node(nodes, "target", row["target_id"], "SPD")
            assay = _add_node(nodes, "assay", row.get("assay_id", f"{row['drug_id']}:{row['target_id']}"), "SPD")
            edges.append(_edge(drug, target, "drug_known_target", "SPD", metadata={"exposure_margin": row.get("exposure_margin", "")}))
            edges.append(_edge(target, assay, "target_assayed", "SPD"))
    nodes_df = pd.DataFrame(nodes.values())
    edges_df = pd.DataFrame(edges).drop_duplicates()
    Path(out_nodes_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_edges_path).parent.mkdir(parents=True, exist_ok=True)
    nodes_df.to_csv(out_nodes_path, index=False)
    edges_df.to_csv(out_edges_path, index=False)
    return nodes_df, edges_df
