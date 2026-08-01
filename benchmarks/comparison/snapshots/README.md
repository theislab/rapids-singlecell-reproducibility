# Result snapshots

Versioned GPU equivalence runs, kept as reviewer-facing evidence. Each directory is named
after the Slurm job that produced it and is never rewritten, so a cited number stays citable.

| Snapshot                                               | Suite               | Node                 | Outcome                       | Status          |
| ------------------------------------------------------ | ------------------- | -------------------- | ----------------------------- | --------------- |
| [`38938766`](38938766)                                 | Expanded, 20 groups | A100-PCIE-40GB       | 16/20 groups, 200/210 metrics | Current         |
| [`../temp-results/38936727`](../temp-results/38936727) | Initial, 15 groups  | A100-PCIE-40GB (MIG) | 14/15 groups, 143/149 metrics | Provenance only |

Only a complete run — every script producing a result record — is promoted here. A partial
rerun marks `execution.json` with `"partial": true` and must not be promoted; use it to iterate
on one method, then take a full run for the snapshot.

Read any snapshot's verdicts against [`../THRESHOLDS.md`](../THRESHOLDS.md), which records the
scientific basis of each criterion and the corrections applied after `38938766`.
