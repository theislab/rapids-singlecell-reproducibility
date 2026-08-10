# Result snapshots

Versioned GPU equivalence runs, kept as reviewer-facing evidence. A snapshot directory is
never rewritten, so a cited number stays citable.

| Snapshot                                     | Suite                                  | Node                 | Outcome                       | Status     |
| -------------------------------------------- | -------------------------------------- | -------------------- | ----------------------------- | ---------- |
| [`2026-08-09-derived`](2026-08-09-derived)   | Comparisons derived from stored arrays | A100-PCIE-40GB (MIG) | 13/20 groups, 165/181 metrics | Current    |
| [`2026-07-31-expanded`](2026-07-31-expanded) | Expanded, 20 groups                    | A100-PCIE-40GB       | 16/20 groups, 200/210 metrics | Superseded |

Only a complete run — every script producing a result record — is promoted here. A partial
rerun marks `execution.json` with `"partial": true` and must not be promoted; use it to iterate
on one method, then take a full run for the snapshot.

Each snapshot's `report/summary.md` lists every failing criterion with the diagnosis behind it.
[`OPEN.md`](../../../OPEN.md) carries the limits of the evidence, and
[`NUMERICAL_VALIDATION.md`](../../../NUMERICAL_VALIDATION.md) the assessment of the `allclose`
criteria against the standard the publication declares.

**Only the current snapshot re-scores.** `evaluate.py` joins records to `criteria.py` by metric
name, and the superseded snapshots predate the systematic naming, so pointing it at them reports
that no criterion applies rather than producing a verdict. Their own recorded numbers, in the table
above and in each directory, remain the citable form.
