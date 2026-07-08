from __future__ import annotations

import importlib
import importlib.util
import json
from pathlib import Path

import pandas as pd


GRAPH_MODEL_BACKENDS = [
    {
        "backend": "pykeen",
        "package": "pykeen",
        "model_family": "knowledge_graph_embedding",
        "example_models": "TransE;DistMult;ComplEx;RotatE;R-GCN via extensions",
        "atlas_role": "Sparse mechanism graph link prediction over drug-target-ADR/pathway edges.",
        "input_requirement": "mechanism_graph_edges.csv with typed source/target nodes and held-out edge split.",
    },
    {
        "backend": "pytorch_geometric",
        "package": "torch_geometric",
        "model_family": "graph_neural_network",
        "example_models": "GCN;GraphSAGE;GAT;R-GCN/HeteroConv",
        "atlas_role": "Heterogeneous graph node/edge representation learning for mechanism support.",
        "input_requirement": "typed node table, edge table, node features, and relation-aware split manifest.",
    },
    {
        "backend": "dgl",
        "package": "dgl",
        "model_family": "graph_neural_network",
        "example_models": "R-GCN;GraphSAGE;GAT;heterograph link prediction",
        "atlas_role": "Alternative heterograph backend for sparse drug-target-ADR mechanism graphs.",
        "input_requirement": "typed node table, edge table, node features, and relation-aware split manifest.",
    },
    {
        "backend": "pyg_hypergraph",
        "package": "torch_geometric",
        "model_family": "hypergraph_neural_network",
        "example_models": "HypergraphConv;heterogeneous hyperedge message passing",
        "atlas_role": "Represent drug-target-ADR, drug-target-pathway-ADR, or drug-target-tissue-ADR tuples as hyperedges instead of flattening them into pair labels.",
        "input_requirement": "hyperedge incidence table, node features, hyperedge labels, and group-aware hyperedge split manifest.",
    },
    {
        "backend": "dgl_hypergraph",
        "package": "dgl",
        "model_family": "hypergraph_neural_network",
        "example_models": "custom DGL heterograph/hypergraph link prediction",
        "atlas_role": "Alternative backend for hyperedge-style ADR mechanism modeling when PyG is not preferred.",
        "input_requirement": "hyperedge incidence table, node features, hyperedge labels, and group-aware hyperedge split manifest.",
    },
    {
        "backend": "kan",
        "package": "torch",
        "package_candidates": "torch;kan;efficient_kan",
        "model_family": "kolmogorov_arnold_network",
        "example_models": "KAN;efficient-KAN;DeepADR-style tabular KAN branch",
        "atlas_role": "Experimental nonlinear tabular expert for sparse ADR mechanism labels; compare against logistic/RF/XGBoost before using in claims.",
        "input_requirement": "same nonleaky tabular expert features plus nested validation/HPO; keep graph evidence out of clean mechanism labels.",
    },
]


def _import_status(candidates: list[str]) -> tuple[list[str], list[str]]:
    available: list[str] = []
    errors: list[str] = []
    for candidate in candidates:
        if importlib.util.find_spec(candidate) is None:
            continue
        try:
            importlib.import_module(candidate)
        except Exception as exc:  # noqa: BLE001 - optional backend import errors must be reported.
            errors.append(f"{candidate}: {type(exc).__name__}: {exc}")
        else:
            available.append(candidate)
    return available, errors


def audit_graph_model_backends() -> pd.DataFrame:
    rows = []
    torch_imports, torch_errors = _import_status(["torch"])
    torch_available = bool(torch_imports)
    for spec in GRAPH_MODEL_BACKENDS:
        package = str(spec["package"])
        candidates = [
            token.strip()
            for token in str(spec.get("package_candidates", package)).split(";")
            if token.strip()
        ]
        present_packages = [candidate for candidate in candidates if importlib.util.find_spec(candidate) is not None]
        available_packages, import_errors = _import_status(candidates)
        available = bool(available_packages) if package != "torch" else torch_available
        rows.append(
            {
                **spec,
                "package_candidates": ";".join(candidates),
                "present_packages": ";".join(present_packages),
                "available_packages": ";".join(available_packages),
                "import_errors": "; ".join(import_errors or torch_errors),
                "available": bool(available),
                "torch_available": bool(torch_available),
                "training_status": "available" if available and (package == "pykeen" or torch_available) else "missing_optional_dependency",
                "production_policy": (
                    "candidate_sparse_graph_model; keep separate from tabular experts until edge-split, relation leakage, and negative-edge sampling are audited"
                ),
            }
        )
    return pd.DataFrame(rows)


def write_graph_model_backend_audit(out_path: str | Path) -> pd.DataFrame:
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame = audit_graph_model_backends()
    frame.to_csv(out, index=False)
    manifest = {
        "out_path": str(out),
        "available_backends": frame.loc[frame["available"], "backend"].astype(str).tolist(),
        "policy": "Graph/KGE backends are optional sparse-mechanism model candidates. They are not used by default tabular expert training.",
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return frame
