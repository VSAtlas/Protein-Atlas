from __future__ import annotations

from pathlib import Path

import pandas as pd


def compute_pair_mechanism_scores(pair_table_path: str | Path, edges_path: str | Path, out_path: str | Path) -> pd.DataFrame:
    pair = pd.read_csv(pair_table_path)
    edges = pd.read_csv(edges_path, low_memory=False)
    edge_set = {(row.source_node_id, row.target_node_id, row.edge_type) for row in edges.itertuples(index=False)}
    by_source: dict[str, set[tuple[str, str]]] = {}
    for row in edges.itertuples(index=False):
        by_source.setdefault(row.source_node_id, set()).add((row.target_node_id, row.edge_type))
    rows: list[dict[str, object]] = []
    for row in pair.itertuples(index=False):
        drug = f"drug:{row.drug_id}"
        target = f"target:{row.target_id}"
        drug_adrs = {dst for dst, edge_type in by_source.get(drug, set()) if edge_type == "drug_has_adr"}
        target_adrs = {dst for dst, edge_type in by_source.get(target, set()) if edge_type == "target_has_safety_evidence"}
        pathways = {dst for dst, edge_type in by_source.get(target, set()) if edge_type == "target_in_pathway"}
        pathway_adrs = {dst2 for pathway in pathways for dst2, edge_type in by_source.get(pathway, set()) if edge_type == "pathway_linked_to_adr"}
        target_pathway_adr_link = bool(pathway_adrs)
        drug_adr_known = bool(drug_adrs)
        target_adr_known = bool(target_adrs)
        drug_target_known = (drug, target, "drug_known_target") in edge_set
        triad_complete = bool(drug_adrs & (target_adrs | pathway_adrs))
        mechanism_path_count = len(drug_adrs & target_adrs) + len(drug_adrs & pathway_adrs)
        score = 1.0 * triad_complete + 0.5 * target_pathway_adr_link + 0.5 * target_adr_known + 0.25 * drug_adr_known
        rows.append(
            {
                "drug_id": row.drug_id,
                "target_id": row.target_id,
                "drug_adr_known": int(drug_adr_known),
                "target_adr_known": int(target_adr_known),
                "target_pathway_adr_link": int(target_pathway_adr_link),
                "drug_target_known": int(drug_target_known),
                "triad_complete": int(triad_complete),
                "mechanism_path_count": mechanism_path_count,
                "mechanism_graph_score": score,
            }
        )
    out = pd.DataFrame(rows)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    return out
