"""Shared Meeko stage review policy."""

from __future__ import annotations

import re
from typing import Literal

MeekoReviewProfile = Literal["pipeline", "publication"]


def _dropped_atom_count(stage: str) -> int:
    match = re.search(r"(?:^|[:_])dropped_atoms=(\d+)", stage)
    if not match:
        return 0
    try:
        return int(match.group(1))
    except ValueError:
        return 0


def meeko_stage_requires_review(
    stage: str, *, profile: MeekoReviewProfile = "pipeline"
) -> bool:
    """Return True when a Meeko-related stage string warrants human review."""

    if profile == "publication":
        return _publication_stage_requires_review(stage)
    return _pipeline_stage_requires_review(stage)


def _publication_stage_requires_review(stage: str) -> bool:
    if not stage:
        return True
    if stage.startswith("internal_oxt_pruned_retry:"):
        return False
    if _dropped_atom_count(stage) > 0:
        return True
    accepted_stages = {
        "modern_clean_input",
        "pdb2pqr_reduce",
        "pdb2pqr_reduce:normalized_noh_retry",
        "pdb2pqr_reduce:normalized_dewatered_noh_retry",
        "normalized_noh_retry",
        "normalized_dewatered_noh_retry",
    }
    if stage in accepted_stages:
        return False
    review_tokens = (
        "fallback",
        "clash_pruned",
        "allow_bad",
        "his_nmap",
        "his_template",
        "legacy",
        "ordered_copy",
        "retry:dropped_atoms",
    )
    return any(token in stage for token in review_tokens)


def _pipeline_stage_requires_review(stage: str) -> bool:
    if stage.startswith("internal_oxt_pruned_retry:"):
        return False
    if _dropped_atom_count(stage) > 0:
        return True
    accepted_stages = {
        "modern_clean_input",
        "modern_his_template_input",
        "his_nmap_retry",
        "normalized_noh_retry",
        "normalized_dewatered_noh_retry",
    }
    if (
        stage in accepted_stages
        or stage.startswith("clash_pruned_safe_retry:")
        or stage.startswith("internal_oxt_pruned_retry:")
    ):
        return False
    review_tokens = (
        "allow_bad",
        "existing",
        "fallback",
        "failed",
        "legacy",
    )
    return any(token in stage for token in review_tokens)
