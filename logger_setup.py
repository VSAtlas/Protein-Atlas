import logging
import os

def setup_logger(pdb_id: str):
    """
    Sets up a logger that logs to 'docked/<pdb_id>/protein.log'.
    Also prints to console.
    """
    log_dir = os.path.join("docked", pdb_id)
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, "protein.log")

    # Remove existing handlers to avoid duplicate logs
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )

    logging.info(f"Logger initialized for PDB: {pdb_id}")
