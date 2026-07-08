from __future__ import annotations

import logging

import post_docking.rescoring.rescoring_scorch as rescoring_scorch


def test_no_combos_after_filters_idempotent_logs_info(caplog) -> None:
    logger = logging.getLogger("test.scorch.idempotent")
    with caplog.at_level(logging.INFO):
        rescoring_scorch._log_no_combos_after_filters(
            logger,
            run_id="rid",
            decoy_prefix="dud",
            pdb_id_filter="TEST",
            variant_filter="HOLO",
            ph_filter="pH7_0",
            idempotent_only=True,
            combos_input=1,
            combos_skipped=1,
        )
    assert any(
        rec.levelno == logging.INFO
        and "no_combos_after_filters_idempotent" in rec.message
        for rec in caplog.records
    )
    assert not any(
        rec.levelno >= logging.WARNING
        and "no_combos_after_filters_idempotent" in rec.message
        for rec in caplog.records
    )


def test_no_combos_after_filters_non_idempotent_logs_warning(caplog) -> None:
    logger = logging.getLogger("test.scorch.non_idempotent")
    with caplog.at_level(logging.INFO):
        rescoring_scorch._log_no_combos_after_filters(
            logger,
            run_id="rid",
            decoy_prefix="dud",
            pdb_id_filter="TEST",
            variant_filter="HOLO",
            ph_filter="pH7_0",
            idempotent_only=False,
            combos_input=0,
            combos_skipped=0,
        )
    assert any(
        rec.levelno == logging.WARNING and "no_combos_after_filters" in rec.message
        for rec in caplog.records
    )
