from __future__ import annotations

import json
import math
import os
from pathlib import Path

HERE = Path(__file__).parent
result_root = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", HERE))
output = Path(os.environ.get("EQUIVALENCE_SUMMARY", HERE / "equivalence.json"))

records = []
invalid_records = []
seen_methods = set()
for path in sorted(result_root.rglob("*.json")):
    if path.resolve() == output.resolve():
        continue
    try:
        record = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        invalid_records.append({"path": str(path), "error": str(error)})
        continue
    if not {"method", "reference_package", "dataset", "tier", "versions", "metrics", "passed"}.issubset(record):
        invalid_records.append({"path": str(path), "error": "Missing required result fields"})
        continue
    if record["method"] in seen_methods:
        invalid_records.append({"path": str(path), "error": f"Duplicate method record: {record['method']}"})
        continue
    if not isinstance(record["metrics"], list) or not record["metrics"]:
        invalid_records.append({"path": str(path), "error": "Metrics must be a non-empty list"})
        continue
    metric_names = set()
    metric_error = None
    for metric in record["metrics"]:
        required = {"metric", "observed", "comparison", "tolerance", "criterion", "passed"}
        if not required.issubset(metric):
            metric_error = "Metric is missing required fields"
            break
        if metric["metric"] in metric_names:
            metric_error = f"Duplicate metric name: {metric['metric']}"
            break
        if metric["comparison"] not in {"<=", ">="}:
            metric_error = f"Unsupported comparison: {metric['comparison']}"
            break
        if not math.isfinite(float(metric["observed"])) or not math.isfinite(float(metric["tolerance"])):
            metric_error = f"Non-finite value in metric: {metric['metric']}"
            break
        metric_names.add(metric["metric"])
    if metric_error:
        invalid_records.append({"path": str(path), "error": metric_error})
        continue
    if bool(record["passed"]) != all(bool(metric["passed"]) for metric in record["metrics"]):
        invalid_records.append({"path": str(path), "error": "Record pass status disagrees with its metrics"})
        continue
    seen_methods.add(record["method"])
    records.append(record)

metrics = [metric for record in records for metric in record["metrics"]]
summary = {
    "passed": bool(records) and not invalid_records and all(record["passed"] for record in records),
    "n_methods": len(records),
    "n_passed_methods": sum(bool(record["passed"]) for record in records),
    "n_metrics": len(metrics),
    "n_passed_metrics": sum(bool(metric["passed"]) for metric in metrics),
    "invalid_records": invalid_records,
    "records": records,
}
output.write_text(json.dumps(summary, indent=2) + "\n")
print(f"Wrote {len(records)} method records and {summary['n_metrics']} metrics to {output}")
if not records:
    raise SystemExit("No result records found. Run comparison scripts first.")
if invalid_records:
    raise SystemExit(f"Could not parse {len(invalid_records)} result record(s).")
if not summary["passed"]:
    raise SystemExit("One or more equivalence thresholds failed.")
