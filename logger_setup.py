import logging, os
from pathlib import Path
from input_and_export_functions import load_inputs, validate_config

def setup_logger(pdb_id: str, cfg: dict | None = None):
    """
    Portable logger.
    Logs to {cfg.get('LOGS_DIR', cfg['DOCKED_DIR'])}/{pdb_id}/protein.log
    Console level honors QUIET_CONSOLE (WARNING if true else INFO).
    """
    if cfg is None:
        cfg = load_inputs(); validate_config(cfg)

    logs_root = Path(cfg.get("LOGS_DIR", cfg["DOCKED_DIR"]))
    log_dir = logs_root / pdb_id
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "protein.log"

    logger = logging.getLogger(f"protein.{pdb_id}")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    fh = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING if bool(cfg.get("QUIET_CONSOLE", False)) else logging.INFO)
    ch.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.propagate = False

    logger.info(f"Logger initialized for PDB: {pdb_id} @ {log_file}")
    return logger
