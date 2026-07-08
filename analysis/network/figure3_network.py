from __future__ import annotations

from pathlib import Path

import pandas as pd


def build_mechanism_network_figure(
    pair_table_path: str | Path,
    mechanism_nodes_path: str | Path,
    mechanism_edges_path: str | Path,
    out_path: str | Path,
    *,
    pbas_pairs_path: str | Path | None = None,
    drug: str | None = None,
    adr: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pair = pd.read_csv(pair_table_path)
    nodes = pd.read_csv(mechanism_nodes_path) if Path(mechanism_nodes_path).exists() else pd.DataFrame()
    edges = pd.read_csv(mechanism_edges_path) if Path(mechanism_edges_path).exists() else pd.DataFrame()
    sub_edges = edges.copy()
    if drug:
        token = str(drug).casefold()
        sub_edges = sub_edges[
            sub_edges.astype(str).apply(lambda col: col.str.casefold().str.contains(token, regex=False)).any(axis=1)
        ]
    if adr:
        token = str(adr).casefold()
        sub_edges = sub_edges[
            sub_edges.astype(str).apply(lambda col: col.str.casefold().str.contains(token, regex=False)).any(axis=1)
        ]
    pair_subset = pair
    if drug and "drug_id" in pair_subset:
        pair_subset = pair_subset[pair_subset["drug_id"].astype(str).str.casefold().str.contains(str(drug).casefold(), regex=False)]
    atlas_edges = _atlas_edges(pair_subset)
    if pbas_pairs_path and Path(pbas_pairs_path).exists():
        atlas_edges = pd.concat([atlas_edges, _pbas_edges(pd.read_csv(pbas_pairs_path), drug)], ignore_index=True)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    case_dir = out.parent.parent / "analysis" / "network_subgraphs" if out.parent.name == "figures" else out.parent
    case_dir.mkdir(parents=True, exist_ok=True)
    case = out.stem.replace("mechanism_network_", "")
    all_edges = pd.concat([sub_edges, atlas_edges], ignore_index=True, sort=False)
    node_ids = set(all_edges.get("source_node_id", pd.Series(dtype=str)).astype(str)) | set(all_edges.get("target_node_id", pd.Series(dtype=str)).astype(str))
    sub_nodes = nodes[nodes.get("node_id", pd.Series(dtype=str)).astype(str).isin(node_ids)] if not nodes.empty else pd.DataFrame({"node_id": sorted(node_ids)})
    sub_nodes.to_csv(case_dir / f"{case}_nodes.csv", index=False)
    all_edges.to_csv(case_dir / f"{case}_edges.csv", index=False)
    _write_graphml(sub_nodes, all_edges, case_dir / f"{case}.graphml")
    _draw_network(sub_nodes, all_edges, out)
    return sub_nodes, all_edges


def _atlas_edges(pair: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in pair.head(200).to_dict("records"):
        drug = f"drug:{row.get('drug_id')}"
        target = f"target:{row.get('target_id')}"
        rows.append(
            {
                "source_node_id": drug,
                "target_node_id": target,
                "edge_type": "Atlas drug-protein prediction",
                "source": "Atlas",
                "atlas_score": row.get("atlas_score", row.get("z_selected", "")),
                "mmgbsa_score": row.get("mmgbsa_score", ""),
                "empirical_p": row.get("combined_empirical_p", row.get("fdr_p_empirical", "")),
                "fdr_q": row.get("fdr_q_value", row.get("fdr_q_target", "")),
                "exposure_plausibility": row.get("exposure_plausibility", ""),
                "tissue_expression": row.get("tissue_expression", ""),
                "target_adr_evidence": row.get("target_adr_evidence", ""),
                "pathway_evidence": row.get("pathway_evidence", ""),
            }
        )
    return pd.DataFrame(rows)


def _pbas_edges(pbas: pd.DataFrame, drug: str | None) -> pd.DataFrame:
    if drug and "drug_id" in pbas:
        pbas = pbas[pbas["drug_id"].astype(str).str.casefold().str.contains(str(drug).casefold(), regex=False)]
    rows = [
        {
            "source_node_id": f"drug:{row.get('drug_id')}",
            "target_node_id": f"target:{row.get('target_id')}",
            "edge_type": "PBAS drug-protein prediction",
            "source": "Sawada_PBAS",
            "pbas_score": row.get("pbas_score", ""),
        }
        for row in pbas.head(200).to_dict("records")
    ]
    return pd.DataFrame(rows)


def _write_graphml(nodes: pd.DataFrame, edges: pd.DataFrame, out_path: Path) -> None:
    lines = ['<?xml version="1.0" encoding="UTF-8"?>', '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">', '<graph edgedefault="directed">']
    for node_id in nodes.get("node_id", pd.Series(dtype=str)).astype(str):
        lines.append(f'<node id="{node_id}"/>')
    for idx, row in enumerate(edges.to_dict("records")):
        lines.append(f'<edge id="e{idx}" source="{row.get("source_node_id", "")}" target="{row.get("target_node_id", "")}"/>')
    lines.extend(["</graph>", "</graphml>"])
    out_path.write_text("\n".join(lines), encoding="utf-8")


def _draw_network(nodes: pd.DataFrame, edges: pd.DataFrame, out_path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5), dpi=160)
    ax.axis("off")
    text = f"Nodes: {len(nodes)}\nEdges: {len(edges)}\nPBAS edges are separated by edge_type/source in exported CSV/GraphML."
    ax.text(0.02, 0.98, text, va="top", ha="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path)
    if out_path.suffix.lower() != ".png":
        fig.savefig(out_path.with_suffix(".png"))
    if out_path.suffix.lower() != ".svg":
        fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)
