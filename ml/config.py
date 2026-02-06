from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


_DEFAULT_CENTER_COLUMNS = ("pocket_center_x", "pocket_center_y", "pocket_center_z")


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
    vina_score: bool = False


@dataclass(frozen=True)
class TrainConfig:
    bigbind_dir: Path
    train_pdb: str = ""
    dataset_split_mode: str = "standard"
    random_seed: int = 42
    max_rows: int | None = None
    splits: tuple[str, ...] = ("train", "val", "test")
    exclude_target_prefixes: tuple[str, ...] = ()
    run_id: str | None = None
    model: ModelConfig = field(default_factory=ModelConfig)
    features: FeaturesConfig = field(default_factory=FeaturesConfig)
    atlas_cfg_path: Path | None = None
    fpocket_center_columns: tuple[str, str, str] = _DEFAULT_CENTER_COLUMNS
    fpocket_variant_column: str | None = "variant"
    fpocket_ph_column: str | None = "pH"
    fpocket_centers_by_pdb: dict[str, tuple[float, float, float]] = field(
        default_factory=dict
    )
    inner_scaffold_folds: int = 5
    source_path: Path | None = None


def _load_json_or_yaml(path: Path) -> dict[str, Any]:
    suffix = path.suffix.lower()
    if suffix in {".json", ".txt"}:
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
        return "standard"
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


def load_config(config_path: str | Path) -> TrainConfig:
    path = Path(config_path).expanduser().resolve()
    payload = _load_json_or_yaml(path)

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
    if isinstance(splits_raw, str):
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

    model = ModelConfig(
        type=str(model_payload.get("type", "logreg")),
        C=float(model_payload.get("C", 1.0)),
        class_weight=class_weight,
        max_iter=int(model_payload.get("max_iter", 1000)),
    )
    features = FeaturesConfig(
        ligand_descriptors=bool(features_payload.get("ligand_descriptors", True)),
        ligand_extra_descriptors=bool(
            features_payload.get("ligand_extra_descriptors", False)
        ),
        ligand_morgan_fp_bits=int(features_payload.get("ligand_morgan_fp_bits", 2048)),
        pocket_fpocket=bool(features_payload.get("pocket_fpocket", True)),
        vina_score=bool(features_payload.get("vina_score", False)),
    )

    center_columns = _normalize_center_columns(payload.get("fpocket_center_columns"))
    variant_column = _coerce_optional_str(
        payload.get("fpocket_variant_column", "variant")
    )
    ph_column = _coerce_optional_str(payload.get("fpocket_ph_column", "pH"))
    centers_by_pdb = _parse_centers_by_pdb(payload.get("fpocket_centers_by_pdb"))
    exclude_target_prefixes = _normalize_prefixes(payload.get("exclude_target_prefixes"))
    inner_scaffold_folds = max(2, int(payload.get("inner_scaffold_folds", 5)))

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
        model=model,
        features=features,
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
    return as_data
