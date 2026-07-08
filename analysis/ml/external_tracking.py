from __future__ import annotations

import os
from pathlib import Path
from typing import Any


def _scalar_metrics(metrics: dict[str, Any]) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, value in metrics.items():
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            continue
        if parsed == parsed:
            out[str(key)] = parsed
    return out


def _simple_params(params: dict[str, Any]) -> dict[str, str | int | float | bool | None]:
    simple: dict[str, str | int | float | bool | None] = {}
    for key, value in params.items():
        if value is None or isinstance(value, str | int | float | bool):
            simple[str(key)] = value
        else:
            simple[str(key)] = str(value)
    return simple


def _infer_ml_root(root: Path) -> Path:
    parts = list(root.resolve().parts)
    if "ml" in parts:
        idx = len(parts) - 1 - parts[::-1].index("ml")
        return Path(*parts[: idx + 1])
    return root


def _infer_run_id(root: Path) -> str:
    parts = list(root.resolve().parts)
    if "data" in parts:
        idx = len(parts) - 1 - parts[::-1].index("data")
        if idx + 1 < len(parts):
            return parts[idx + 1]
    return os.environ.get("ATLAS_RUN_ID", "atlas")


def log_external_trackers(
    *,
    model_dir: str | Path,
    params: dict[str, Any],
    metrics: dict[str, Any],
    dataset_version_id: str | None,
    claim_status: str | None,
) -> dict[str, dict[str, Any]]:
    """Best-effort logging for Atlas local registry plus optional third-party trackers.

    MLflow is enabled by default with a local run-scoped file store. Set
    ATLAS_DISABLE_MLFLOW=1 or pass the wrapper --no-mlflow flag to suppress it.
    """
    root = Path(model_dir)
    status: dict[str, dict[str, Any]] = {}
    mlflow_enabled = os.environ.get("ATLAS_DISABLE_MLFLOW") != "1"
    if mlflow_enabled:
        status["mlflow"] = _log_mlflow(root, params, metrics, dataset_version_id, claim_status)
    if os.environ.get("ATLAS_WANDB_PROJECT"):
        status["wandb"] = _log_wandb(root, params, metrics, dataset_version_id, claim_status)
    if os.environ.get("ATLAS_CLEARML_PROJECT"):
        status["clearml"] = _log_clearml(root, params, metrics, dataset_version_id, claim_status)
    return status


def _log_mlflow(
    root: Path,
    params: dict[str, Any],
    metrics: dict[str, Any],
    dataset_version_id: str | None,
    claim_status: str | None,
) -> dict[str, Any]:
    try:
        import mlflow
    except Exception as exc:
        return {"status": "skipped", "reason": f"mlflow import failed: {exc}"}
    try:
        tracking_uri = os.environ.get("ATLAS_MLFLOW_TRACKING_URI")
        if not tracking_uri:
            tracking_uri = f"file:{_infer_ml_root(root) / 'mlruns'}"
        mlflow.set_tracking_uri(tracking_uri)
        registry_uri = os.environ.get("ATLAS_MLFLOW_REGISTRY_URI")
        if registry_uri:
            mlflow.set_registry_uri(registry_uri)
        mlflow.set_experiment(os.environ.get("ATLAS_MLFLOW_EXPERIMENT", f"atlas/{_infer_run_id(root)}/ml"))
        with mlflow.start_run(run_name=root.name):
            mlflow.log_params(_simple_params(params))
            mlflow.log_metrics(_scalar_metrics(metrics))
            mlflow.set_tags(
                {
                    "dataset_version_id": dataset_version_id or "",
                    "claim_status": claim_status or "",
                    "atlas_model_dir": str(root),
                    "atlas_ml_root": str(_infer_ml_root(root)),
                }
            )
            mlflow.log_artifacts(str(root))
            registered_model: dict[str, Any] | None = None
            if os.environ.get("ATLAS_MLFLOW_REGISTER_MODEL") == "1":
                try:
                    import pickle

                    with (root / "trained_model.pkl").open("rb") as handle:
                        model = pickle.load(handle)
                    model_name = os.environ.get("ATLAS_MLFLOW_REGISTERED_MODEL_NAME") or f"atlas-{_infer_run_id(root)}-{root.name}"
                    model_info = mlflow.sklearn.log_model(
                        sk_model=model,
                        artifact_path="registered_model",
                        registered_model_name=model_name,
                    )
                    registered_model = {
                        "status": "registered",
                        "name": model_name,
                        "model_uri": getattr(model_info, "model_uri", None),
                        "registry_uri": registry_uri or "default",
                    }
                except Exception as exc:  # pragma: no cover - optional registry integration
                    registered_model = {"status": "failed", "reason": str(exc)}
            run_id = mlflow.active_run().info.run_id if mlflow.active_run() else None
        result = {"status": "logged", "run_id": run_id, "tracking_uri": tracking_uri}
        if registry_uri:
            result["registry_uri"] = registry_uri
        if registered_model is not None:
            result["registered_model"] = registered_model
        return result
    except Exception as exc:
        return {"status": "failed", "reason": str(exc)}


def _log_wandb(
    root: Path,
    params: dict[str, Any],
    metrics: dict[str, Any],
    dataset_version_id: str | None,
    claim_status: str | None,
) -> dict[str, Any]:
    try:
        import wandb
    except Exception as exc:
        return {"status": "skipped", "reason": f"wandb import failed: {exc}"}
    try:
        mode = os.environ.get("WANDB_MODE") or ("online" if os.environ.get("ATLAS_WANDB_ONLINE") == "1" else "offline")
        run = wandb.init(
            project=os.environ["ATLAS_WANDB_PROJECT"],
            name=root.name,
            dir=str(root),
            mode=mode,
            config={**_simple_params(params), "dataset_version_id": dataset_version_id, "claim_status": claim_status},
            reinit=True,
        )
        wandb.log(_scalar_metrics(metrics))
        run.finish()
        return {"status": "logged", "mode": mode, "run_id": getattr(run, "id", None)}
    except Exception as exc:
        return {"status": "failed", "reason": str(exc)}


def _log_clearml(
    root: Path,
    params: dict[str, Any],
    metrics: dict[str, Any],
    dataset_version_id: str | None,
    claim_status: str | None,
) -> dict[str, Any]:
    try:
        from clearml import Task
    except Exception as exc:
        return {"status": "skipped", "reason": f"clearml import failed: {exc}"}
    try:
        task = Task.init(
            project_name=os.environ["ATLAS_CLEARML_PROJECT"],
            task_name=os.environ.get("ATLAS_CLEARML_TASK", root.name),
            output_uri=False,
        )
        task.connect({**_simple_params(params), "dataset_version_id": dataset_version_id, "claim_status": claim_status})
        logger = task.get_logger()
        for key, value in _scalar_metrics(metrics).items():
            logger.report_scalar(title="metrics", series=key, value=value, iteration=0)
        task.close()
        return {"status": "logged", "task_id": task.id}
    except Exception as exc:
        return {"status": "failed", "reason": str(exc)}
