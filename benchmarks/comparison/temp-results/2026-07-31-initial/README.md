# Temporary GPU equivalence results

> Superseded by [`../../snapshots/2026-07-31-expanded`](../../snapshots/2026-07-31-expanded), the
> first complete run of
> the expanded 20-group suite. This directory is kept only as provenance for the earlier
> 15-group suite; see [`../../snapshots/README.md`](../../snapshots/README.md) for the index.

This directory contains a temporary reference snapshot from a run on 2026-07-31.
It covers the initial 15 structured method groups and predates the core Scanpy and
end-to-end biological additions. It should not be interpreted as a result from the expanded
reviewer suite.

## Outcome

- Methods passing: 14/15
- Metrics passing: 143/149
- Runtime: 8 minutes 23 seconds
- GPU: NVIDIA A100-PCIE-40GB, 20 GB MIG partition
- Driver: 575.57.08
- CUDA: 12.9
- rapids-singlecell: 0.16.1
- scanpy: 1.12.3
- squidpy: 1.8.3
- pertpy: 1.1.1
- decoupler: 2.2.0

| Comparison          | Result | Metrics |
| ------------------- | ------ | ------: |
| bbknn_scrublet      | PASS   |     5/5 |
| calculate_niche     | FAIL   |     3/9 |
| clustering_extended | PASS   |     4/4 |
| co_occurrence       | PASS   |     4/4 |
| decoupler_methods   | PASS   |   18/18 |
| distance            | PASS   |   18/18 |
| embeddings_extended | PASS   |     6/6 |
| guide_assignment    | PASS   |     3/3 |
| ingest_cell_cycle   | PASS   |     6/6 |
| ligrec              | PASS   |     4/4 |
| mixscale            | PASS   |     2/2 |
| mixscape            | PASS   |     5/5 |
| rank_genes_groups   | PASS   |   60/60 |
| spatial_autocorr    | PASS   |     4/4 |
| sqrt                | PASS   |     1/1 |

The CellCharter portion of `calculate_niche` passed. Neighborhood and UTAG clustering did not meet the configured equivalence thresholds; their exact measurements are recorded in [equivalence.json](equivalence.json).

## Files

- [equivalence.json](equivalence.json) contains every metric, tolerance, and pass/fail decision.
- [execution.json](execution.json) records the exit status of every comparison script.
- [environment.txt](environment.txt) records the pinned Python environment.

This snapshot is not updated automatically. Automatic reporting requires GPU-backed CI on a self-hosted GPU runner. That CI should run the suite, retain the JSON and environment files as artifacts, and publish the aggregate result as a job or pull-request summary.
