from __future__ import annotations

import logging
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import prep_ligands.prep_ligands_common as ligprep_common
import prep_ligands.prep_ligands_reporting as ligprep_reporting
from config.tool_resolver import resolve_tool
from config.runtime_config import load_config, validate_config
from path_router import make_paths
from prep_ligands.prep_ligands_common import get_short_path_name
from prep_ligands.prep_ligands_microstates import (
    load_microstate_registry,
    save_microstate_registry,
)

logger = logging.getLogger(__name__)


@dataclass(slots=True)
class PrepToolchain:
    mgltools_python_short: str
    prepare_script_short: str
    obabel_exe_short: str


def resolve_prep_toolchain(cfg: Dict[str, Any]) -> PrepToolchain:
    obabel_exe = str(
        resolve_tool(cfg, "OPENBABEL_PATH", "obabel").get("resolved_path", "") or ""
    ).strip()
    if obabel_exe and not os.environ.get("BABEL_DATADIR"):
        data_dir = Path(obabel_exe).resolve().parent / "data"
        if data_dir.exists():
            os.environ["BABEL_DATADIR"] = str(data_dir)

    return PrepToolchain(
        mgltools_python_short="",
        prepare_script_short="",
        obabel_exe_short=get_short_path_name(obabel_exe) if obabel_exe else "",
    )


def init_prep_workspace(
    output_dir: Path,
    *,
    source_link_target: Optional[Path] = None,
    create_reference_dirs: bool = False,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    if create_reference_dirs:
        for dirname in ("intermediates", "reference", "quarantine"):
            (output_dir / dirname).mkdir(parents=True, exist_ok=True)
    if source_link_target is not None:
        link = output_dir / "intermediates_src"
        try:
            if link.exists() or link.is_symlink():
                if link.is_dir() and not link.is_symlink():
                    shutil.rmtree(link)
                else:
                    link.unlink()
            link.symlink_to(source_link_target.resolve(), target_is_directory=True)
            logging.info("[intermediates] symlinked source -> intermediates_src")
        except Exception as exc:
            logging.debug("intermediates link skipped: %s", exc)
    return output_dir / "ligand_prep_status.tsv"


@dataclass(slots=True)
class MicrostateState:
    enabled: bool
    library_name: str
    library_out_dir: Path
    directory: Optional[Path] = None
    registry: Optional[dict] = None
    index: Optional[dict] = None
    alias_index: Optional[Dict[Tuple[str, str, float], dict]] = None
    dirty: bool = False

    def save_if_needed(self, *, force: bool = False) -> None:
        if not self.enabled or self.registry is None:
            return
        microstate_list = self.registry.get("microstates") or []
        alias_count = sum(len(entry.get("aliases") or []) for entry in microstate_list)
        if microstate_list:
            logger.info(
                "prep_ligands: microstate_dedup summary library=%s out_dir=%s microstates=%d aliases=%d",
                self.library_name,
                self.library_out_dir,
                len(microstate_list),
                alias_count,
            )
        else:
            logger.warning(
                "prep_ligands: microstate_dedup summary library=%s out_dir=%s microstates=%d aliases=%d (EMPTY REGISTRY)",
                self.library_name,
                self.library_out_dir,
                len(microstate_list),
                alias_count,
            )
        if self.dirty or force:
            save_microstate_registry(self.library_out_dir, self.registry)
            self.dirty = False

    def get(self, key: str, default: Any = None) -> Any:
        return {
            "microstates_dir": self.directory,
            "microstate_registry": self.registry,
            "microstate_index": self.index,
            "alias_index": self.alias_index,
            "microstate_registry_dirty": self.dirty,
        }.get(key, default)

    def __setitem__(self, key: str, value: Any) -> None:
        if key == "microstate_registry_dirty":
            self.dirty = bool(value)
            return
        raise KeyError(key)


@dataclass(slots=True)
class BulkContext:
    cfg: Dict[str, Any]
    force: bool
    library_name: str
    library_hint: str
    ligand_extracted_dir: Path
    ligands_mol2_dir: Path
    output_ligands_dir: Path
    prepped_ligands_dir: Path
    status_log: Path
    mgltools_python_short: str
    prepare_script_short: str
    obabel_exe_short: str
    in_sdf_env: Optional[str]
    in_pdb_dir_env: Optional[str]
    rename_prefix: str
    rename_pad: int
    rename_start: int
    rename_force: bool
    rename_active: bool
    is_fda_library: bool
    microstate_dedup: bool
    use_ph_subdirs: bool
    ph_values: Optional[List[float]]
    eff_ph_values: List[float]
    use_rdkit_for_3d: bool
    obabel_threads: int
    obabel_timeout_s: int
    chunk_size: int
    ligprep_ph: float
    keep_nonpolar_h: int
    max_heavy_atoms: int
    min_atoms_for_docking: int
    min_parent_heavy: int
    min_tors_dof: int
    microstates: MicrostateState = field(repr=False)

    def save_microstates(self, *, force: bool = False) -> None:
        self.microstates.save_if_needed(force=force)

    def _compat_map(self) -> Dict[str, Any]:
        return {
            "cfg": self.cfg,
            "force": self.force,
            "library_name": self.library_name,
            "library_hint": self.library_hint,
            "ligand_extracted_dir": self.ligand_extracted_dir,
            "ligands_mol2_dir": self.ligands_mol2_dir,
            "output_ligands_dir": self.output_ligands_dir,
            "prepped_ligands_dir": self.prepped_ligands_dir,
            "library_out_dir": self.prepped_ligands_dir,
            "status_log": self.status_log,
            "mgltools_python_short": self.mgltools_python_short,
            "prepare_script_short": self.prepare_script_short,
            "obabel_exe_short": self.obabel_exe_short,
            "in_sdf_env": self.in_sdf_env,
            "in_pdb_dir_env": self.in_pdb_dir_env,
            "rename_prefix": self.rename_prefix,
            "rename_pad": self.rename_pad,
            "rename_start": self.rename_start,
            "rename_force": self.rename_force,
            "rename_active": self.rename_active,
            "is_fda_library": self.is_fda_library,
            "microstate_dedup": self.microstate_dedup,
            "use_ph_subdirs": self.use_ph_subdirs,
            "ph_values": self.ph_values,
            "eff_ph_values": self.eff_ph_values,
            "USE_RDKIT_FOR_3D": self.use_rdkit_for_3d,
            "OBABEL_THREADS": self.obabel_threads,
            "OBABEL_TIMEOUT_S": self.obabel_timeout_s,
            "CHUNK_SIZE": self.chunk_size,
            "LIGPREP_PH": self.ligprep_ph,
            "KEEP_NONPOLAR_H": self.keep_nonpolar_h,
            "MAX_HEAVY_ATOMS": self.max_heavy_atoms,
            "MIN_ATOMS_FOR_DOCKING": self.min_atoms_for_docking,
            "MIN_PARENT_HEAVY": self.min_parent_heavy,
            "MIN_TORS_DOF": self.min_tors_dof,
            "microstates_dir": self.microstates.directory,
            "microstate_registry": self.microstates.registry,
            "microstate_index": self.microstates.index,
            "alias_index": self.microstates.alias_index,
            "microstate_registry_dirty": self.microstates.dirty,
            "maybe_save_microstate_registry": self.save_microstates,
        }

    def __getitem__(self, key: str) -> Any:
        compat = self._compat_map()
        if key not in compat:
            raise KeyError(key)
        return compat[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self._compat_map().get(key, default)

    def __setitem__(self, key: str, value: Any) -> None:
        attr_map = {
            "ligands_mol2_dir": "ligands_mol2_dir",
            "microstate_registry_dirty": None,
        }
        if key == "microstate_registry_dirty":
            self.microstates.dirty = bool(value)
            return
        attr = attr_map.get(key)
        if attr is None:
            raise KeyError(key)
        setattr(self, attr, value)


@dataclass(slots=True)
class BulkInitArgs:
    force: bool
    microstate_dedup: bool
    ph_values: Optional[List[float]]
    in_sdf: Optional[Path]
    in_sdf_dir: Optional[Path]
    in_pdb_dir: Optional[Path]
    mol2_dir: Optional[Path]
    out_pdbqt_dir: Optional[Path]
    status_log: Optional[Path]
    root_dir: Optional[Path]


def _env_bool(name: str, default: Any) -> bool:
    return str(os.environ.get(name, default)).strip().lower() in {"1", "true", "yes", "y"}


def _env_int(name: str, default: Any) -> int:
    return int(str(os.environ.get(name, default)).strip())


def _env_float(name: str, default: Any) -> float:
    return float(str(os.environ.get(name, default)).strip())


def sanitize_ligand_name_for_filename(name: str) -> str:
    import re

    name = name.strip().replace(" ", "_")
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", name)[:80] or "ligand"


def relative_to_output(path: Path, output_dir: Path) -> str:
    try:
        return str(path.relative_to(output_dir))
    except Exception:
        return path.name


def ph_label(ph: float) -> str:
    return f"pH{ph:.1f}".replace(".", "_")


def init_bulk_context(args: BulkInitArgs) -> BulkContext:
    force = args.force or _env_bool("LIGPREP_FORCE", False)
    print("Starting ligand preparation")
    if args.microstate_dedup:
        logger.info(
            "prep_ligands: microstate_dedup=True (stage 4: microstate registry + canonical PDBQT shadow)"
        )

    cfg = load_config()
    validate_config(cfg)
    pdb_token = (
        os.environ.get("PDB_ID")
        or cfg.get("PDB_ID")
        or cfg.get("TARGET_PDB")
        or cfg.get("PDB")
        or cfg.get("INPUT_PDB")
        or cfg.get("PDB_FILE")
        or ""
    )
    pdb_id = Path(str(pdb_token)).stem.upper() if pdb_token else "LIGPREP"
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    env_values = {
        key: (os.environ.get(env_name, "") or "").strip() or None
        for key, env_name in {
            "in_sdf_env": "LIGPREP_IN_SDF",
            "in_sdf_dir_env": "LIGPREP_IN_SDF_DIR",
            "in_pdb_dir_env": "LIGPREP_IN_PDB_DIR",
            "mol2_dir_env": "LIGPREP_MOL2_DIR",
            "out_dir_env": "LIGPREP_OUT_DIR",
            "status_log_env": "LIGPREP_STATUS_LOG",
        }.items()
    }
    explicit_paths = {
        "in_sdf_env": args.in_sdf,
        "in_sdf_dir_env": args.in_sdf_dir,
        "in_pdb_dir_env": args.in_pdb_dir,
        "mol2_dir_env": args.mol2_dir,
        "out_dir_env": args.out_pdbqt_dir,
        "status_log_env": args.status_log,
    }
    for key, value in explicit_paths.items():
        if value is not None:
            env_values[key] = str(Path(value).expanduser().resolve())

    rename_prefix = (os.environ.get("LIGPREP_RENAME_PREFIX", "") or "").strip()
    try:
        rename_pad = max(1, int(os.environ.get("LIGPREP_RENAME_PAD", "5")))
    except Exception:
        rename_pad = 5
    try:
        rename_start = max(1, int(os.environ.get("LIGPREP_RENAME_START", "1")))
    except Exception:
        rename_start = 1
    rename_force = _env_bool("LIGPREP_RENAME_FORCE", False)
    rename_active = bool(rename_prefix or rename_force)

    output_ligands_dir = (
        Path(env_values["out_dir_env"]).resolve()
        if env_values["out_dir_env"]
        else paths.prepped_ligands_dir
    )
    prepped_ligands_dir = output_ligands_dir
    library_hint = output_ligands_dir.name.lower()
    if args.root_dir is not None:
        output_ligands_dir = args.root_dir.resolve()
        prepped_ligands_dir = output_ligands_dir
        library_hint = prepped_ligands_dir.name.lower()

    library_env = (os.environ.get("LIGPREP_LIBRARY", "") or "").strip().lower()
    library_base = library_env or library_hint or "ligprep"
    if not library_env and env_values["in_sdf_env"]:
        try:
            in_sdf_path = Path(env_values["in_sdf_env"]).resolve()
            parts = list(in_sdf_path.parts)
            extracted_idx = parts.index("extracted_ligands") if "extracted_ligands" in parts else -1
            detected = parts[extracted_idx + 1] if extracted_idx >= 0 and extracted_idx + 1 < len(parts) else ""
            library_base = (detected or in_sdf_path.stem).lower()
        except Exception:
            pass
    if not library_env and args.root_dir is None:
        os.environ["LIGPREP_LIBRARY"] = library_base

    ligand_extracted_dir = paths.ligand_output_dir
    if args.root_dir is not None and not env_values["in_sdf_dir_env"]:
        project_root = next((parent.parent for parent in args.root_dir.parents if parent.name == "prepped_ligands"), None)
        project_root = project_root or args.root_dir.parent.parent
        candidate_extracted = project_root / "extracted_ligands" / library_base
        if candidate_extracted.is_dir():
            ligand_extracted_dir = candidate_extracted
            logger.info("prep_ligands: using extracted_ligands source path=%s", candidate_extracted)
    elif env_values["in_sdf_dir_env"]:
        ligand_extracted_dir = Path(env_values["in_sdf_dir_env"]).resolve()

    ligands_mol2_dir = (
        Path(env_values["mol2_dir_env"]).resolve()
        if env_values["mol2_dir_env"]
        else Path(paths.ligands_mol2_dir) / library_base
    )
    if "prepped_ligands" in str(ligands_mol2_dir):
        suggested = Path(str(ligands_mol2_dir).replace("prepped_ligands", "ligands_mol2")).resolve()
        logging.warning("[compat] LIGANDS_MOL2_DIR points at prepped_ligands; redirecting to %s", suggested)
        ligands_mol2_dir = suggested

    ligprep_common.MALFORMED_DIR = prepped_ligands_dir
    ligprep_reporting.MALFORMED_DIR = prepped_ligands_dir
    init_prep_workspace(output_ligands_dir)
    ligands_mol2_dir.mkdir(parents=True, exist_ok=True)

    microstates_dir = prepped_ligands_dir / "microstates"
    if args.microstate_dedup:
        microstates_dir.mkdir(parents=True, exist_ok=True)
        registry, index = load_microstate_registry(prepped_ligands_dir, library_base)
        alias_index: Dict[Tuple[str, str, float], dict] = {}
        for entry in registry.get("microstates", []) or []:
            for alias in entry.get("aliases", []) or []:
                try:
                    alias_index[(str(alias.get("ligand_stem", "")), str(alias.get("ph_label", "")), float(alias.get("ph_value", 0.0)))] = entry
                except Exception:
                    continue
        microstates = MicrostateState(
            enabled=True,
            library_name=library_base,
            library_out_dir=prepped_ligands_dir,
            directory=microstates_dir,
            registry=registry,
            index=index,
            alias_index=alias_index,
        )
    else:
        microstates = MicrostateState(enabled=False, library_name=library_base, library_out_dir=prepped_ligands_dir)

    use_rdkit_for_3d = _env_bool("USE_RDKIT_FOR_3D", cfg.get("USE_RDKIT_FOR_3D", True))
    obabel_threads = _env_int("OBABEL_THREADS", cfg.get("OBABEL_THREADS", 50))
    obabel_timeout_s = _env_int("OBABEL_TIMEOUT_S", cfg.get("OBABEL_TIMEOUT_S", 900))
    chunk_size = _env_int("LIGPREP_CHUNK_SIZE", cfg.get("LIGPREP_CHUNK_SIZE", 200))
    ligprep_ph = _env_float("LIGPREP_PH", cfg.get("LIGPREP_PH", 7.4))
    keep_nonpolar_h = _env_int("KEEP_NONPOLAR_H", cfg.get("KEEP_NONPOLAR_H", 1))
    max_heavy_atoms = _env_int("MAX_HEAVY_ATOMS", cfg.get("MAX_HEAVY_ATOMS", 1200))
    min_atoms_for_docking = _env_int("MIN_ATOMS_FOR_DOCKING", cfg.get("MIN_ATOMS_FOR_DOCKING", 5))
    min_parent_heavy = _env_int("MIN_PARENT_HEAVY", cfg.get("MIN_PARENT_HEAVY", 8))
    min_tors_dof = _env_int("MIN_TORS_DOF", cfg.get("MIN_TORS_DOF", 0))

    eff_ph_values = [float(x) for x in args.ph_values] if args.ph_values else [ligprep_ph]
    ph_source = "python" if args.ph_values is not None else "config"
    print(f"[ligprep] ph_schedule source={ph_source} values={','.join(f'{ph:.1f}' for ph in eff_ph_values)} multiple={len(eff_ph_values) > 1}")
    logger.info("prep_ligands: pH schedule eff_ph_values=%s", eff_ph_values)

    toolchain = resolve_prep_toolchain(cfg)
    for label, value in {
        "EXTRACTED_LIGANDS_DIR": ligand_extracted_dir,
        "LIGANDS_MOL2_DIR": ligands_mol2_dir,
        "PREPPED_LIGANDS_DIR": output_ligands_dir,
    }.items():
        if not str(value).strip():
            raise RuntimeError(
                f"Ligand library path check failed: {label} is missing/empty "
                f"(value='{value}'). Set it in config.txt and run `atlas --doctor`."
            )
    if not ligand_extracted_dir.exists():
        raise FileNotFoundError(
            "Ligand library input check failed: EXTRACTED_LIGANDS_DIR does not exist at "
            f"{ligand_extracted_dir}. Generate/provide ligand inputs or update config.txt, "
            "then run `atlas --doctor`."
        )

    status_log_basename = env_values["status_log_env"] or cfg.get(
        "LIGAND_STATUS_LOG_BASENAME", "ligand_prep_status.tsv"
    )
    status_log_path = (
        Path(env_values["status_log_env"]).resolve()
        if env_values["status_log_env"] and Path(env_values["status_log_env"]).suffix
        else output_ligands_dir / status_log_basename
    )

    print(
        "[paths.effective]"
        f" in_sdf={env_values['in_sdf_env'] or 'None'}"
        f" in_sdf_dir={ligand_extracted_dir}"
        f" in_pdb_dir={env_values['in_pdb_dir_env'] or 'None'}"
        f" mol2_dir={ligands_mol2_dir}"
        f" out_pdbqt_dir={output_ligands_dir}"
        f" status_log={(env_values['status_log_env'] or cfg.get('LIGAND_STATUS_LOG_BASENAME', 'ligand_prep_status.tsv'))}"
    )
    print(f"[ligprep] sdf_dir={ligands_mol2_dir / '_rdkit_embedded_sdf'} mol2_dir={ligands_mol2_dir} pdbqt_dir={output_ligands_dir}")

    return BulkContext(
        cfg=cfg,
        force=force,
        library_name=library_base,
        library_hint=library_hint,
        ligand_extracted_dir=ligand_extracted_dir,
        ligands_mol2_dir=ligands_mol2_dir,
        output_ligands_dir=output_ligands_dir,
        prepped_ligands_dir=prepped_ligands_dir,
        status_log=status_log_path,
        mgltools_python_short=toolchain.mgltools_python_short,
        prepare_script_short=toolchain.prepare_script_short,
        obabel_exe_short=toolchain.obabel_exe_short,
        in_sdf_env=env_values["in_sdf_env"],
        in_pdb_dir_env=env_values["in_pdb_dir_env"],
        rename_prefix=rename_prefix,
        rename_pad=rename_pad,
        rename_start=rename_start,
        rename_force=rename_force,
        rename_active=rename_active,
        is_fda_library=library_hint == "fda",
        microstate_dedup=args.microstate_dedup,
        use_ph_subdirs=args.ph_values is not None,
        ph_values=args.ph_values,
        eff_ph_values=eff_ph_values,
        use_rdkit_for_3d=use_rdkit_for_3d,
        obabel_threads=obabel_threads,
        obabel_timeout_s=obabel_timeout_s,
        chunk_size=chunk_size,
        ligprep_ph=ligprep_ph,
        keep_nonpolar_h=keep_nonpolar_h,
        max_heavy_atoms=max_heavy_atoms,
        min_atoms_for_docking=min_atoms_for_docking,
        min_parent_heavy=min_parent_heavy,
        min_tors_dof=min_tors_dof,
        microstates=microstates,
    )
