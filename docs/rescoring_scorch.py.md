# rescoring_scorch.py changes

## Summary

1. Tightened chunking policy under queue pressure (smaller chunks, stronger packing).
2. Lowered tail split trigger and runtime rechunk stale threshold for earlier tail balancing.
3. Updated elastic thread grants to favor more concurrent low-thread jobs under fragmentation.

## Why

Utilization benchmarks showed fragmentation and dispatch gaps dominate idle core time. These policy changes bias toward better fit and less stranded tail capacity.

## Expected effect

1. Higher chunk concurrency when cores are fragmented.
2. Faster tail drain with earlier rechunking.
3. More stable throughput under mixed-size SCORCH task loads.
