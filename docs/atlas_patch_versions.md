# Atlas Patch Versions

Patch labels use `AtlasvMAJOR.MINOR.PATCH` followed by a short description. They identify repository patches and analysis interface changes; they are not claims that datasets or scientific results are immutable releases.

## Atlasv0.0.01 - Context-preserving Phase 1 PK ingestion

Adds a one-command PK refresh for AtlasSPD Phase 1. The workflow preserves every source-specific PK context row, writes a separately selected representative record, joins external values under `pk_context_*` columns, and does not recompute `spd_exposure_label`. SPD is primary; NCATS FRDB and PK-DB availability is recorded explicitly; DailyMed/openFDA requests are serial, cached, and low-confidence until numeric extraction is reviewed; DrugBank remains optional BYOL.
