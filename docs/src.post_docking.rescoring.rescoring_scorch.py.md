# Notes: `src/post_docking/rescoring/rescoring_scorch.py`

- Added `TASK_*` bench event normalization through `_emit_task_event(...)` so SCORCH tasks consistently emit `want_cores`, `min_cores`, and `granted_cores`.
- Added adaptive chunking for SCORCH ligand batches (`_adaptive_chunk_size`) using runtime scheduler free-core/queue snapshots.
- Added tail splitting for large late-stage chunks (`_tail_split_allowed_chunks` + tail queue split pass) to reduce end-of-run stranded cores.
- Added `est_duration_sec` propagation into scheduler admission calls for SCORCH materialize/score/post tasks to improve duration-aware packing.
- Replaced bulk `ThreadPoolExecutor` submission with a managed pending/running admission loop so runtime rechunking can occur while work is in flight.
- Added tail runtime rechunk trigger when free cores remain stranded and pending chunks are uniformly large.
- Added capped tail hedging for long-running straggler chunks, including loser cleanup and hedge-output promotion events.
- Added per-job thread elasticity (`_elastic_scorch_threads`) so SCORCH scoring scales down under fragmentation and scales up when spare cores are abundant.
- Hardened SCORCH event lifecycle so `TASK_END` is emitted even when materialization/scoring raises exceptions; error details are now included in the event payload.
