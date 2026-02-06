from __future__ import annotations

import json
import os
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


_DEFAULT_CENTER_COLUMNS = ("pocket_center_x", "pocket_center_y", "pocket_center_z")
_ENV_VAR_PATTERN = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")

_TXT_KEY_ALIASES: dict[str, str] = {
    "BIGBIND_DIR": "bigbind_dir",
    "TRAIN_PDB": "train_pdb",
    "DATASET_SPLIT_MODE": "dataset_split_mode",
    "SPLITS": "splits",
    "EXCLUDE_TARGET_PREFIXES": "exclude_target_prefixes",
    "RANDOM_SEED": "random_seed",
    "MAX_ROWS": "max_rows",
    "RUN_ID": "run_id",
    "ATLAS_CFG_PATH": "atlas_cfg_path",
    "FPOCKET_CENTER_COLUMNS": "fpocket_center_columns",
    "FPOCKET_VARIANT_COLUMN": "fpocket_variant_column",
    "FPOCKET_PH_COLUMN": "fpocket_ph_column",
    "INNER_SCAFFOLD_FOLDS": "inner_scaffold_folds",
    "MODEL_TYPE": "model.type",
    "MODEL_C": "model.c",
    "MODEL_CLASS_WEIGHT": "model.class_weight",
    "MODEL_MAX_ITER": "model.max_iter",
    "MODEL_FAMILY": "model_family",
    "MODEL_PARAMS": "model_params",
    "FEATURES_LIGAND_DESCRIPTORS": "features.ligand_descriptors",
    "FEATURES_LIGAND_EXTRA_DESCRIPTORS": "features.ligand_extra_descriptors",
    "FEATURES_LIGAND_MORGAN_FP_BITS": "features.ligand_morgan_fp_bits",
    "FEATURES_POCKET_FPOCKET": "features.pocket_fpocket",
    "FEATURES_POCKET_FEATURES": "features.pocket_features",
    "FEATURES_VINA_SCORE": "features.vina_score",
    "CALIBRATION_ENABLED": "calibration.enabled",
    "CALIBRATION_METHOD": "calibration.method",
    "CALIBRATION_CV_FOLDS": "calibration.cv_folds",
    "CALIBRATION_SEED": "calibration.seed",
    "METRICS_REPORT_BRIER": "metrics.report_brier",
    "METRICS_REPORT_ECE": "metrics.report_ece",
    "METRICS_ECE_BINS": "metrics.ece_bins",
    "REGISTRY_ENABLED": "registry.enabled",
}

_TXT_LIST_KEYS = {
    "splits",
    "exclude_target_prefixes",
    "fpocket_center_columns",
    "automl.feature_variants",
    "automl.model_families",
}


@dataclass(frozen=True)
class ModelConfig:
    type: str = "logreg"
    C: float = 1.0
    class_weight: str | None = "balanced"
    max_iter: int = 1000


@dataclass(frozen=True)
class FeaturesConfig:
    ligand_descriptors: bool = True
    ligand_extra_descriptors: bool = False
    ligand_morgan_fp_bits: int = 2048
    pocket_fpocket: bool = True
    pocket_features: bool = False
    vina_score: bool = False


@dataclass(frozen=True)
class RankConfig:
    enabled: bool = False
    group_key: str = "target_pocket"


@dataclass(frozen=True)
class HardNegativesConfig:
    enabled: bool = False
    policy: str = "property"
    per_active_k: int = 5
    within_group: bool = True
    max_candidates: int = 2000
    ratio: float = 1.0
    seed: int = 42


@dataclass(frozen=True)
class LabelsConfig:
    smoothing_eps: float = 0.0
    use_sample_weights: bool = False
    weight_cap: float = 5.0


@dataclass(frozen=True)
class ProteinFamilySplitConfig:
    identity_threshold: float = 0.3
    seed: int = 42


@dataclass(frozen=True)
class PocketSimilaritySplitConfig:
    n_clusters: int = 5
    seed: int = 42


@dataclass(frozen=True)
class StressSplitConfig:
    mode: str = "standard"
    protein_family: ProteinFamilySplitConfig = field(default_factory=ProteinFamilySplitConfig)
    pocket_similarity: PocketSimilaritySplitConfig = field(
        default_factory=PocketSimilaritySplitConfig
    )


@dataclass(frozen=True)
class PipelineStageConfig:
    model_family: str = "xgboost"
    model_params: dict[str, Any] = field(default_factory=dict)
    keep_top_pct: float | None = None
    keep_prob_ge: float | None = None
    keep_within_group: bool = True
    group_key: str = "target_pocket"
    calibration_enabled: bool = False
    calibration_method: str = "sigmoid"
    calibration_cv_folds: int = 5
    calibration_seed: int = 42


@dataclass(frozen=True)
class PipelineConfig:
    enabled: bool = False
    stage1: PipelineStageConfig = field(default_factory=PipelineStageConfig)
    stage2: PipelineStageConfig = field(
        default_factory=lambda: PipelineStageConfig(
            model_family="lightgbm_rank",
            group_key="target_pocket",
        )
    )


@dataclass(frozen=True)
class DiagnosticsConfig:
    adversarial_validation: bool = True
    drop_feature_tests: bool = True
    permutation_importance: bool = True
    permutation_max_rows: int = 5000


@dataclass(frozen=True)
class TrainConfig:
    bigbind_dir: Path
    train_pdb: str = ""
    dataset_split_mode: str = "legacy"
    random_seed: int = 42
    max_rows: int | None = None
    splits: tuple[str, ...] = ("train", "val", "test")
    exclude_target_prefixes: tuple[str, ...] = ()
    run_id: str | None = None
    model_family: str = "logreg"
    model_params: dict[str, Any] = field(default_factory=dict)
    model: ModelConfig = field(default_factory=ModelConfig)
    features: FeaturesConfig = field(default_factory=FeaturesConfig)
    rank: RankConfig = field(default_factory=RankConfig)
    hard_negatives: HardNegativesConfig = field(default_factory=HardNegativesConfig)
    labels: LabelsConfig = field(default_factory=LabelsConfig)
    stress_splits: StressSplitConfig = field(default_factory=StressSplitConfig)
    pipeline: PipelineConfig = field(default_factory=PipelineConfig)
    diagnostics: DiagnosticsConfig = field(default_factory=DiagnosticsConfig)
    calibration_enabled: bool = False
    calibration_method: str = "sigmoid"
    calibration_cv_folds: int = 5
    calibration_seed: int = 42
    metrics_report_brier: bool = True
    metrics_report_ece: bool = True
    metrics_ece_bins: int = 15
    registry_enabled: bool = True
    atlas_cfg_path: Path | None = None
    fpocket_center_columns: tuple[str, str, str] = _DEFAULT_CENTER_COLUMNS
    fpocket_variant_column: str | None = "variant"
    fpocket_ph_column: str | None = "pH"
    fpocket_centers_by_pdb: dict[str, tuple[float, float, float]] = field(
        default_factory=dict
    )
    inner_scaffold_folds: int = 5
    source_path: Path | None = None


def _strip_inline_comment(text: str) -> str:
    in_single = False
    in_double = False
    for idx, char in enumerate(text):
        if char == "'" and not in_double:
            in_single = not in_single
            continue
        if char == '"' and not in_single:
            in_double = not in_double
            continue
        if char == "#" and not in_single and not in_double:
            return text[:idx].rstrip()
    return text.strip()


def _expand_template_vars(value: str, variables: dict[str, str]) -> str:
    def _replace(match: re.Match[str]) -> str:
        key = match.group(1)
        if key in variables:
            return variables[key]
        return os.environ.get(key, match.group(0))

    expanded = value
    for _ in range(10):
        updated = _ENV_VAR_PATTERN.sub(_replace, expanded)
        if updated == expanded:
            break
        expanded = updated
    return expanded


def _coerce_text_value(raw: str, *, key_path: str) -> Any:
    lowered = raw.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    if lowered in {"null", "none"}:
        return None

    if key_path in _TXT_LIST_KEYS and raw and not raw.lstrip().startswith("["):
        return [part.strip() for part in raw.split(",") if part.strip()]

    if raw and raw[0] in "[{":
        try:
            return json.loads(raw)
        except Exception:
            pass

    if re.fullmatch(r"[+-]?\d+", raw):
        try:
            return int(raw)
        except Exception:
            pass
    if re.fullmatch(r"[+-]?(?:\d+\.\d*|\d*\.\d+)(?:[eE][+-]?\d+)?", raw):
        try:
            return float(raw)
        except Exception:
            pass

    return raw


def _canonical_txt_key(raw_key: str) -> str:
    stripped = raw_key.strip()
    if not stripped:
        return stripped
    alias = _TXT_KEY_ALIASES.get(stripped)
    if alias is not None:
        return alias
    alias = _TXT_KEY_ALIASES.get(stripped.upper())
    if alias is not None:
        return alias
    if "." in stripped:
        return ".".join(part.strip().lower() for part in stripped.split(".") if part.strip())
    return stripped.lower()


def _set_nested(payload: dict[str, Any], dotted_key: str, value: Any) -> None:
    if "." not in dotted_key:
        payload[dotted_key] = value
        return
    parts = [part for part in dotted_key.split(".") if part]
    if not parts:
        return
    cursor: dict[str, Any] = payload
    for part in parts[:-1]:
        existing = cursor.get(part)
        if not isinstance(existing, dict):
            existing = {}
            cursor[part] = existing
        cursor = existing
    cursor[parts[-1]] = value


def _load_key_value_payload(path: Path) -> dict[str, Any]:
    raw_text = path.read_text(encoding="utf-8")
    stripped = raw_text.lstrip()
    if stripped.startswith("{"):
        loaded = json.loads(raw_text)
        return loaded if isinstance(loaded, dict) else {}

    payload: dict[str, Any] = {}
    variables: dict[str, str] = {}
    for line in raw_text.splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        if "=" not in text:
            continue
        raw_key, raw_value = text.split("=", 1)
        key = _canonical_txt_key(raw_key)
        value_text = _strip_inline_comment(raw_value.strip())
        if not value_text:
            parsed_value: Any = ""
        else:
            if (
                len(value_text) >= 2
                and value_text[0] == value_text[-1]
                and value_text[0] in {"'", '"'}
            ):
                value_text = value_text[1:-1]
            expanded = _expand_template_vars(value_text, variables)
            parsed_value = _coerce_text_value(expanded, key_path=key)

        _set_nested(payload, key, parsed_value)

        as_text = "" if parsed_value is None else str(parsed_value)
        variables[raw_key.strip()] = as_text
        variables[raw_key.strip().upper()] = as_text
        variables[key] = as_text
        variables[key.upper()] = as_text
    return payload


def _load_json_or_yaml(path: Path) -> dict[str, Any]:
    suffix = path.suffix.lower()
    if suffix == ".txt":
        return _load_key_value_payload(path)
    if suffix == ".json":
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    if suffix in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore[import-untyped]
        except Exception as exc:  # pragma: no cover - optional dependency
            raise RuntimeError(
                "YAML config requested but PyYAML is not installed. "
                "Use JSON config or install PyYAML."
            ) from exc
        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
        return loaded or {}
    raise ValueError(f"Unsupported config extension: {path.suffix}")


def _coerce_optional_int(value: Any) -> int | None:
    if value in (None, "", "none", "null"):
        return None
    return int(value)


def _coerce_optional_str(value: Any, *, default: str | None = None) -> str | None:
    if value is None:
        return default
    text = str(value).strip()
    if text.lower() in {"", "none", "null"}:
        return None
    return text


def _resolve_path(raw_path: str, *, base_dir: Path) -> Path:
    candidate = Path(raw_path).expanduser()
    if candidate.is_absolute():
        return candidate
    return (base_dir / candidate).resolve()


def _normalize_center_columns(raw: Any) -> tuple[str, str, str]:
    if raw is None:
        return _DEFAULT_CENTER_COLUMNS
    if not isinstance(raw, (list, tuple)) or len(raw) != 3:
        raise ValueError(
            "Config key 'fpocket_center_columns' must be a list/tuple with 3 values."
        )
    columns = tuple(str(col).strip() for col in raw)
    if not all(columns):
        raise ValueError("Config key 'fpocket_center_columns' cannot contain blanks.")
    return columns  # type: ignore[return-value]


def _parse_center_triplet(raw: Any) -> tuple[float, float, float] | None:
    if raw is None:
        return None
    if not isinstance(raw, (list, tuple)) or len(raw) != 3:
        return None
    try:
        return (float(raw[0]), float(raw[1]), float(raw[2]))
    except Exception:
        return None


def _parse_centers_by_pdb(raw: Any) -> dict[str, tuple[float, float, float]]:
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise ValueError("Config key 'fpocket_centers_by_pdb' must be a dict.")
    parsed: dict[str, tuple[float, float, float]] = {}
    for key, value in raw.items():
        pdb_id = str(key).strip().lower()
        if not pdb_id:
            continue
        center = _parse_center_triplet(value)
        if center is not None:
            parsed[pdb_id] = center
    return parsed


def _normalize_prefixes(raw: Any) -> tuple[str, ...]:
    if raw is None:
        return ()
    if isinstance(raw, str):
        values = [part.strip() for part in raw.split(",")]
    elif isinstance(raw, (list, tuple)):
        values = [str(item).strip() for item in raw]
    else:
        raise ValueError(
            "Config key 'exclude_target_prefixes' must be a string or list/tuple of strings."
        )
    normalized = tuple(value for value in values if value)
    return normalized


def _normalize_dataset_split_mode(raw: Any) -> str:
    if raw is None:
        return "legacy"
    mode = str(raw).strip().lower()
    if mode in {"standard", "train_val_test"}:
        return "standard"
    if mode in {"legacy", "train_pdb_holdout"}:
        return "legacy"
    raise ValueError(
        "Config key 'dataset_split_mode' must be one of: "
        "standard, train_val_test, legacy, train_pdb_holdout."
    )


def _load_default_atlas_cfg() -> dict[str, Any]:
    try:
        from analysis.dud_eval_orchestrate import _load_default_cfg

        cfg = _load_default_cfg() or {}
        if isinstance(cfg, dict):
            return cfg
    except Exception:
        pass

    try:
        from input_and_export_functions import load_config as load_atlas_config

        cfg = load_atlas_config("config.txt") or {}
        return cfg if isinstance(cfg, dict) else {}
    except Exception:
        return {}


def _load_atlas_cfg_from_path(path: Path) -> dict[str, Any]:
    suffix = path.suffix.lower()
    if suffix in {".json", ".yaml", ".yml"}:
        loaded = _load_json_or_yaml(path)
        if not isinstance(loaded, dict):
            raise ValueError(f"atlas_cfg_path must load to a dict: {path}")
        merged = _load_default_atlas_cfg()
        merged.update(loaded)
        return merged

    try:
        from input_and_export_functions import load_config as load_atlas_config

        loaded = load_atlas_config(config_path=str(path), base_dir=path.parent) or {}
        return loaded if isinstance(loaded, dict) else {}
    except TypeError:
        loaded = load_atlas_config(str(path)) or {}
        return loaded if isinstance(loaded, dict) else {}


def load_atlas_cfg(config: TrainConfig) -> dict[str, Any]:
    if config.atlas_cfg_path is None:
        return _load_default_atlas_cfg()
    return _load_atlas_cfg_from_path(config.atlas_cfg_path)


def load_config_payload(config_path: str | Path) -> dict[str, Any]:
    path = Path(config_path).expanduser().resolve()
    payload = _load_json_or_yaml(path)
    if not isinstance(payload, dict):
        raise ValueError(f"Config file must contain a dict/object: {path}")
    return payload


def load_config(config_path: str | Path) -> TrainConfig:
    path = Path(config_path).expanduser().resolve()
    payload = load_config_payload(path)

    if "bigbind_dir" not in payload:
        raise KeyError("Missing required config key: bigbind_dir")
    dataset_split_mode = _normalize_dataset_split_mode(payload.get("dataset_split_mode"))
    train_pdb_value = str(payload.get("train_pdb", "")).strip()
    if dataset_split_mode == "legacy" and not train_pdb_value:
        raise KeyError(
            "Missing required config key: train_pdb when dataset_split_mode is legacy."
        )

    model_payload = payload.get("model", {}) or {}
    features_payload = payload.get("features", {}) or {}
    model_params_payload = payload.get("model_params", {}) or {}
    calibration_payload = payload.get("calibration", {}) or {}
    metrics_payload = payload.get("metrics", {}) or {}
    registry_payload = payload.get("registry", {}) or {}
    rank_payload = payload.get("rank", {}) or {}
    hard_negative_payload = payload.get("hard_negatives", {}) or {}
    labels_payload = payload.get("labels", {}) or {}
    splits_payload = payload.get("stress_splits", payload.get("split_regime", {})) or {}
    pipeline_payload = payload.get("pipeline", {}) or {}
    diagnostics_payload = payload.get("diagnostics", {}) or {}

    class_weight_raw = model_payload.get("class_weight", "balanced")
    if isinstance(class_weight_raw, str) and class_weight_raw.strip().lower() in {
        "",
        "none",
        "null",
    }:
        class_weight = None
    else:
        class_weight = class_weight_raw

    splits_raw = payload.get("splits", ("train", "val", "test"))
    if isinstance(splits_raw, dict):
        split_names_raw = splits_raw.get("names", ("train", "val", "test"))
        if isinstance(split_names_raw, str):
            splits = tuple(
                part.strip() for part in split_names_raw.split(",") if part.strip()
            )
        else:
            splits = tuple(str(part).strip() for part in split_names_raw if str(part).strip())
        if not splits_payload:
            splits_payload = splits_raw
    elif isinstance(splits_raw, str):
        splits = tuple(part.strip() for part in splits_raw.split(",") if part.strip())
    else:
        splits = tuple(str(part).strip() for part in splits_raw if str(part).strip())
    if not splits:
        raise ValueError("Config key 'splits' cannot be empty.")

    atlas_cfg_path_raw = payload.get("atlas_cfg_path")
    atlas_cfg_path = (
        None
        if atlas_cfg_path_raw in (None, "", "none", "null")
        else _resolve_path(str(atlas_cfg_path_raw), base_dir=path.parent)
    )

    model_c_raw = model_payload.get("C", model_payload.get("c", 1.0))
    model = ModelConfig(
        type=str(model_payload.get("type", "logreg")),
        C=float(model_c_raw),
        class_weight=class_weight,
        max_iter=int(model_payload.get("max_iter", 1000)),
    )

    model_family_raw = payload.get("model_family", model.type)
    model_family = str(model_family_raw or "logreg").strip().lower()
    if model_family in {"lgbm"}:
        model_family = "lightgbm"
    if model_family in {"xgb"}:
        model_family = "xgboost"

    if isinstance(model_params_payload, dict):
        model_params = {str(k): v for k, v in model_params_payload.items()}
    else:
        model_params = {}
    if model_family in {"logreg", "logistic", "logistic_regression"}:
        model_params.setdefault("C", model.C)
        model_params.setdefault("class_weight", model.class_weight)
        model_params.setdefault("max_iter", model.max_iter)
    for legacy_key in ("penalty", "solver"):
        if legacy_key in model_payload and legacy_key not in model_params:
            model_params[legacy_key] = model_payload[legacy_key]

    features = FeaturesConfig(
        ligand_descriptors=bool(features_payload.get("ligand_descriptors", True)),
        ligand_extra_descriptors=bool(
            features_payload.get("ligand_extra_descriptors", False)
        ),
        ligand_morgan_fp_bits=int(features_payload.get("ligand_morgan_fp_bits", 2048)),
        pocket_fpocket=bool(features_payload.get("pocket_fpocket", True)),
        pocket_features=bool(features_payload.get("pocket_features", False)),
        vina_score=bool(features_payload.get("vina_score", False)),
    )
    rank_cfg = RankConfig(
        enabled=bool(rank_payload.get("enabled", False)),
        group_key=str(rank_payload.get("group_key", "target_pocket")).strip().lower(),
    )
    hard_negatives_cfg = HardNegativesConfig(
        enabled=bool(hard_negative_payload.get("enabled", False)),
        policy=str(hard_negative_payload.get("policy", "property")).strip().lower(),
        per_active_k=max(1, int(hard_negative_payload.get("per_active_k", 5))),
        within_group=bool(hard_negative_payload.get("within_group", True)),
        max_candidates=max(1, int(hard_negative_payload.get("max_candidates", 2000))),
        ratio=max(0.0, float(hard_negative_payload.get("ratio", 1.0))),
        seed=int(hard_negative_payload.get("seed", payload.get("random_seed", 42))),
    )
    labels_cfg = LabelsConfig(
        smoothing_eps=float(labels_payload.get("smoothing_eps", 0.0)),
        use_sample_weights=bool(labels_payload.get("use_sample_weights", False)),
        weight_cap=max(1.0, float(labels_payload.get("weight_cap", 5.0))),
    )
    protein_family_cfg = ProteinFamilySplitConfig(
        identity_threshold=float(
            (splits_payload.get("protein_family", {}) or {}).get("identity_threshold", 0.3)
        ),
        seed=int((splits_payload.get("protein_family", {}) or {}).get("seed", payload.get("random_seed", 42))),
    )
    pocket_similarity_cfg = PocketSimilaritySplitConfig(
        n_clusters=max(2, int((splits_payload.get("pocket_similarity", {}) or {}).get("n_clusters", 5))),
        seed=int((splits_payload.get("pocket_similarity", {}) or {}).get("seed", payload.get("random_seed", 42))),
    )
    stress_split_cfg = StressSplitConfig(
        mode=str(splits_payload.get("mode", "standard")).strip().lower(),
        protein_family=protein_family_cfg,
        pocket_similarity=pocket_similarity_cfg,
    )

    stage1_payload = pipeline_payload.get("stage1", {}) or {}
    stage2_payload = pipeline_payload.get("stage2", {}) or {}
    pipeline_cfg = PipelineConfig(
        enabled=bool(pipeline_payload.get("enabled", False)),
        stage1=PipelineStageConfig(
            model_family=str(stage1_payload.get("model_family", "xgboost")).strip().lower(),
            model_params=dict(stage1_payload.get("model_params", {}) or {}),
            keep_top_pct=(
                None
                if stage1_payload.get("keep_top_pct", None) in (None, "", "none", "null")
                else float(stage1_payload.get("keep_top_pct"))
            ),
            keep_prob_ge=(
                None
                if stage1_payload.get("keep_prob_ge", None) in (None, "", "none", "null")
                else float(stage1_payload.get("keep_prob_ge"))
            ),
            keep_within_group=bool(stage1_payload.get("keep_within_group", True)),
            group_key=str(stage1_payload.get("group_key", "target_pocket")).strip().lower(),
            calibration_enabled=bool(stage1_payload.get("calibration_enabled", False)),
            calibration_method=str(stage1_payload.get("calibration_method", "sigmoid")).strip().lower(),
            calibration_cv_folds=max(2, int(stage1_payload.get("calibration_cv_folds", 5))),
            calibration_seed=int(stage1_payload.get("calibration_seed", payload.get("random_seed", 42))),
        ),
        stage2=PipelineStageConfig(
            model_family=str(stage2_payload.get("model_family", "lightgbm_rank")).strip().lower(),
            model_params=dict(stage2_payload.get("model_params", {}) or {}),
            keep_top_pct=None,
            keep_prob_ge=None,
            keep_within_group=bool(stage2_payload.get("keep_within_group", True)),
            group_key=str(stage2_payload.get("group_key", "target_pocket")).strip().lower(),
            calibration_enabled=bool(stage2_payload.get("calibration_enabled", False)),
            calibration_method=str(stage2_payload.get("calibration_method", "sigmoid")).strip().lower(),
            calibration_cv_folds=max(2, int(stage2_payload.get("calibration_cv_folds", 5))),
            calibration_seed=int(stage2_payload.get("calibration_seed", payload.get("random_seed", 42))),
        ),
    )
    diagnostics_cfg = DiagnosticsConfig(
        adversarial_validation=bool(diagnostics_payload.get("adversarial_validation", True)),
        drop_feature_tests=bool(diagnostics_payload.get("drop_feature_tests", True)),
        permutation_importance=bool(diagnostics_payload.get("permutation_importance", True)),
        permutation_max_rows=max(100, int(diagnostics_payload.get("permutation_max_rows", 5000))),
    )

    center_columns = _normalize_center_columns(payload.get("fpocket_center_columns"))
    variant_column = _coerce_optional_str(
        payload.get("fpocket_variant_column", "variant")
    )
    ph_column = _coerce_optional_str(payload.get("fpocket_ph_column", "pH"))
    centers_by_pdb = _parse_centers_by_pdb(payload.get("fpocket_centers_by_pdb"))
    exclude_target_prefixes = _normalize_prefixes(payload.get("exclude_target_prefixes"))
    inner_scaffold_folds = max(2, int(payload.get("inner_scaffold_folds", 5)))
    calibration_enabled = bool(calibration_payload.get("enabled", False))
    calibration_method = str(calibration_payload.get("method", "sigmoid")).strip().lower()
    if calibration_method not in {"sigmoid", "isotonic"}:
        raise ValueError("Config calibration.method must be 'sigmoid' or 'isotonic'.")
    calibration_cv_folds = max(2, int(calibration_payload.get("cv_folds", 5)))
    calibration_seed = int(calibration_payload.get("seed", payload.get("random_seed", 42)))

    metrics_report_brier = bool(metrics_payload.get("report_brier", True))
    metrics_report_ece = bool(metrics_payload.get("report_ece", True))
    metrics_ece_bins = max(2, int(metrics_payload.get("ece_bins", 15)))
    registry_enabled = bool(registry_payload.get("enabled", True))

    return TrainConfig(
        bigbind_dir=_resolve_path(str(payload["bigbind_dir"]), base_dir=path.parent),
        train_pdb=train_pdb_value,
        dataset_split_mode=dataset_split_mode,
        random_seed=int(payload.get("random_seed", 42)),
        max_rows=_coerce_optional_int(payload.get("max_rows")),
        splits=splits,
        exclude_target_prefixes=exclude_target_prefixes,
        run_id=(
            str(payload["run_id"]).strip()
            if payload.get("run_id") not in (None, "")
            else None
        ),
        model_family=model_family,
        model_params=model_params,
        model=model,
        features=features,
        rank=rank_cfg,
        hard_negatives=hard_negatives_cfg,
        labels=labels_cfg,
        stress_splits=stress_split_cfg,
        pipeline=pipeline_cfg,
        diagnostics=diagnostics_cfg,
        calibration_enabled=calibration_enabled,
        calibration_method=calibration_method,
        calibration_cv_folds=calibration_cv_folds,
        calibration_seed=calibration_seed,
        metrics_report_brier=metrics_report_brier,
        metrics_report_ece=metrics_report_ece,
        metrics_ece_bins=metrics_ece_bins,
        registry_enabled=registry_enabled,
        atlas_cfg_path=atlas_cfg_path,
        fpocket_center_columns=center_columns,
        fpocket_variant_column=variant_column,
        fpocket_ph_column=ph_column,
        fpocket_centers_by_pdb=centers_by_pdb,
        inner_scaffold_folds=inner_scaffold_folds,
        source_path=path,
    )


def config_to_dict(config: TrainConfig) -> dict[str, Any]:
    as_data = asdict(config)
    as_data["bigbind_dir"] = str(config.bigbind_dir)
    as_data["atlas_cfg_path"] = (
        str(config.atlas_cfg_path) if config.atlas_cfg_path else None
    )
    as_data["fpocket_center_columns"] = list(config.fpocket_center_columns)
    as_data["splits"] = list(config.splits)
    as_data["exclude_target_prefixes"] = list(config.exclude_target_prefixes)
    as_data["source_path"] = str(config.source_path) if config.source_path else None
    as_data["calibration"] = {
        "enabled": bool(config.calibration_enabled),
        "method": str(config.calibration_method),
        "cv_folds": int(config.calibration_cv_folds),
        "seed": int(config.calibration_seed),
    }
    as_data["metrics"] = {
        "report_brier": bool(config.metrics_report_brier),
        "report_ece": bool(config.metrics_report_ece),
        "ece_bins": int(config.metrics_ece_bins),
    }
    as_data["registry"] = {"enabled": bool(config.registry_enabled)}
    as_data["rank"] = {
        "enabled": bool(config.rank.enabled),
        "group_key": str(config.rank.group_key),
    }
    as_data["hard_negatives"] = {
        "enabled": bool(config.hard_negatives.enabled),
        "policy": str(config.hard_negatives.policy),
        "per_active_k": int(config.hard_negatives.per_active_k),
        "within_group": bool(config.hard_negatives.within_group),
        "max_candidates": int(config.hard_negatives.max_candidates),
        "ratio": float(config.hard_negatives.ratio),
        "seed": int(config.hard_negatives.seed),
    }
    as_data["labels"] = {
        "smoothing_eps": float(config.labels.smoothing_eps),
        "use_sample_weights": bool(config.labels.use_sample_weights),
        "weight_cap": float(config.labels.weight_cap),
    }
    as_data["stress_splits"] = {
        "mode": str(config.stress_splits.mode),
        "protein_family": {
            "identity_threshold": float(config.stress_splits.protein_family.identity_threshold),
            "seed": int(config.stress_splits.protein_family.seed),
        },
        "pocket_similarity": {
            "n_clusters": int(config.stress_splits.pocket_similarity.n_clusters),
            "seed": int(config.stress_splits.pocket_similarity.seed),
        },
    }
    as_data["pipeline"] = {
        "enabled": bool(config.pipeline.enabled),
        "stage1": asdict(config.pipeline.stage1),
        "stage2": asdict(config.pipeline.stage2),
    }
    as_data["diagnostics"] = {
        "adversarial_validation": bool(config.diagnostics.adversarial_validation),
        "drop_feature_tests": bool(config.diagnostics.drop_feature_tests),
        "permutation_importance": bool(config.diagnostics.permutation_importance),
        "permutation_max_rows": int(config.diagnostics.permutation_max_rows),
    }
    return as_data
