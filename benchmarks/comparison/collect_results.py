from __future__ import annotations

import json
import os
from pathlib import Path


HERE = Path(__file__).parent
result_root = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", HERE))
output = Path(os.environ.get("EQUIVALENCE_SUMMARY", HERE / "equivalence.json"))

records = []
for path in sorted(result_root.rglob("*.json")):
    if path.resolve() == output.resolve():
        continue
    record = json.loads(path.read_text())
    if {"method", "metrics", "passed"}.issubset(record):
        records.append(record)

summary = {
    "passed": bool(records) and all(record["passed"] for record in records),
    "n_methods": len(records),
    "n_metrics": sum(len(record["metrics"]) for record in records),
    "records": records,
}
output.write_text(json.dumps(summary, indent=2) + "\n")
print(f"Wrote {len(records)} method records and {summary['n_metrics']} metrics to {output}")
if not records:
    raise SystemExit("No result records found. Run comparison scripts first.")
if not summary["passed"]:
    raise SystemExit("One or more equivalence thresholds failed.")
