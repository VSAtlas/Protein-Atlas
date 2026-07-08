from __future__ import annotations

import logging
from pathlib import Path

from protein_prep.pdb_fixer_runtime import (
    fix_element_columns_in_file,
    rules_version,
    scan_helium_counts,
    scan_helium_counts_with_hits,
)


def _post_write_element_guard(step: str, pdb_path: str | Path) -> None:
    """
    Normalize element columns after each receptor write step and log helium status.
    """
    p = Path(pdb_path)
    try:
        fixed = fix_element_columns_in_file(p, dst_path=p, rewrite_atoms=False)
    except Exception as e:
        logging.warning(
            "[helium] stage=post_write step=%s file=%s He->H=? residual_He=? note=elemfix_error:%s",
            step,
            p,
            e,
        )
        return

    txt = p.read_text(encoding="utf-8", errors="ignore")
    he_count = scan_helium_counts(txt)
    logging.info(
        "[helium] stage=post_write step=%s file=%s He->H=%s residual_He=%d",
        step,
        p,
        fixed if isinstance(fixed, int) else -1,
        he_count,
    )


def _meeko_preflight_or_fail(pdb_input: str | Path, work_dir: str | Path) -> Path:
    """
    Ensure the exact PDB Meeko will read contains no 'He'.
    Saves copy to work/meeko_input_pre_sanitize.pdb for forensics.
    """
    src = Path(pdb_input).resolve()
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)

    try:
        fix_element_columns_in_file(src, dst_path=src, rewrite_atoms=False)
    except Exception as e:
        logging.warning(
            "[helium] stage=preflight file=%s note=elemfix_exception:%s", src, e
        )

    txt = src.read_text(encoding="utf-8", errors="ignore")
    he_count, first_hits = scan_helium_counts_with_hits(txt, max_hits=5)

    # Save what Meeko will actually see
    prefile = work / "meeko_input_pre_sanitize.pdb"
    prefile.write_text(txt, encoding="utf-8")

    if he_count > 0:
        logging.error(
            "[helium] stage=preflight file=%s He_count=%d rules=%s",
            src,
            he_count,
            rules_version(),
        )
        for hit in first_hits:
            logging.error("[helium] offender %s", hit)
        raise RuntimeError("helium_preflight_failed")

    return src
