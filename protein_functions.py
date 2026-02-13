import logging
import os
import subprocess

import activesite

logger = logging.getLogger(__name__)


def detect_active_site(cleaned_pdb):
    exists = os.path.exists(cleaned_pdb)
    size = os.path.getsize(cleaned_pdb) if exists else -1
    logger.info(
        "[detect_active_site] entry cleaned_pdb=%s exists=%s size=%s",
        cleaned_pdb,
        exists,
        size,
    )
    label = "unknown"
    base = os.path.basename(cleaned_pdb)
    if "nolig" in base.lower():
        label = "nolig"
    elif "cleaned" in base.lower():
        label = "cleaned"
    logger.debug("[detect_active_site] cleaned_pdb_label=%s", label)

    center = box_size = source = None
    try:
        center, box_size, source = activesite.main(cleaned_pdb)
    except Exception as exc:
        logger.warning(
            "[detect_active_site] activesite.main failed pdb=%s err=%s",
            cleaned_pdb,
            exc,
        )
        return None, None, None

    logger.info(
        "[detect_active_site] exit center=%s box=%s source=%s cleaned_pdb=%s",
        center,
        box_size,
        source,
        cleaned_pdb,
    )
    if center is None:
        print("❌ Active site detection failed.")
    return center, box_size, source


def convert_to_pdbqt(cleaned_pdb, receptor_pdbqt, mgltools_python, prepare_script):
    try:
        subprocess.run(
            [mgltools_python, prepare_script, "-r", cleaned_pdb, "-o", receptor_pdbqt],
            check=True,
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to convert protein: {e}")
        return False
