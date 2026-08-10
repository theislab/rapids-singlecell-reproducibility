# Result snapshots

Versioned GPU equivalence runs, kept as reviewer-facing evidence. A snapshot directory is
never rewritten, so a cited number stays citable.

| Snapshot                                   | Suite                                  | Node                 | Outcome                       | Status  |
| ------------------------------------------ | -------------------------------------- | -------------------- | ----------------------------- | ------- |
| [`2026-08-09-derived`](2026-08-09-derived) | Comparisons derived from stored arrays | A100-PCIE-40GB (MIG) | 15/20 groups, 149/161 metrics | Current |

Only a complete run — every script producing a result record — is promoted here. A partial
rerun marks `execution.json` with `"partial": true` and must not be promoted; use it to iterate
on one method, then take a full run for the snapshot.

Each snapshot's `report/summary.md` lists every failing criterion with the diagnosis behind it.
[`OPEN.md`](../../../OPEN.md) carries the limits of the evidence, and
[`NUMERICAL_VALIDATION.md`](../../../NUMERICAL_VALIDATION.md) the assessment of the `allclose`
criteria against the standard the publication declares.

Promoting a snapshot means committing its `report/logs/`, which the repository-wide `*.log`
rule otherwise excludes: add a matching negation in `.gitignore`, once the logs are confirmed
to carry no host paths. Only container runs qualify, since their paths are all `/repro`.
