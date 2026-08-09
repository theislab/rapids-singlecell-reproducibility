"""Evaluate stored measurements against their criteria. No GPU, no rerun.

The comparison scripts measure; this decides. They are separate programs on purpose:
measuring costs a GPU and ten minutes, whereas deciding costs milliseconds, so any
question of the form "what would the verdict be if..." is answered here against records
that already exist.

Criteria come from `criteria.py`, never from the record. Whatever a script wrote about
thresholds is a measurement-time artefact and is overwritten here, so an old result
directory — including a committed snapshot — can be re-scored under today's criteria.

    python benchmarks/comparison/evaluate.py
    python benchmarks/comparison/evaluate.py --results snapshots/<name>/results
    python benchmarks/comparison/evaluate.py --criteria what-if.toml

The optional TOML file overrides thresholds without touching `criteria.py`. Keys are
`"<method>.<metric>"`, because a metric name alone is ambiguous — the same name carries
different thresholds in different method groups:

    ["biological_pipeline_pbmc3k.umap.cross_embedding_knn_overlap"]
    tolerance = 0.4
    basis = "what-if: at the reseeded-CPU baseline"
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from criteria import criterion_for, diagnosis_for

HERE = Path(__file__).parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--results",
        type=Path,
        default=Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", HERE / "results")),
        help="Directory of per-method measurement records.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path(os.environ.get("EQUIVALENCE_SUMMARY", HERE / "equivalence.json")),
        help="Where to write the aggregated evaluation.",
    )
    parser.add_argument("--criteria", type=Path, help='TOML file overriding thresholds, keyed "<method>.<metric>".')
    parser.add_argument("--execution", type=Path, help="Execution record; enables report rendering.")
    parser.add_argument("--report-dir", type=Path, help="Render the reviewer report here.")
    return parser.parse_args()


def load_overrides(path: Path | None) -> dict:
    if path is None:
        return {}
    overrides = tomllib.loads(path.read_text())
    for metric, entry in overrides.items():
        if "tolerance" not in entry:
            raise SystemExit(f"Criteria override for {metric!r} must set a tolerance")
    return overrides


def decide(metric: dict) -> bool:
    """Apply the metric's criterion. The single place a verdict is formed.

    A non-finite measurement fails. NaN is not close to anything, and an unrepresentable
    result is not evidence of agreement — treating it as a bad record instead would throw
    away the whole method group over one metric.
    """
    if not math.isfinite(float(metric["observed"])):
        return False
    if metric["comparison"] == "<=":
        return float(metric["observed"]) <= float(metric["tolerance"])
    return float(metric["observed"]) >= float(metric["tolerance"])


def validate(record: dict) -> str | None:
    """Reject a malformed measurement record rather than silently scoring it.

    Only measurements are required here. Anything a script may still write about
    thresholds is ignored: criteria.py is the authority on what a criterion is.
    """
    if not {"method", "reference_package", "dataset", "tier", "versions", "metrics"}.issubset(record):
        return "Missing required result fields"
    if not isinstance(record["metrics"], list) or not record["metrics"]:
        return "Metrics must be a non-empty list"
    seen = set()
    for metric in record["metrics"]:
        if not {"metric", "observed"}.issubset(metric):
            return "Metric is missing required fields"
        if metric["metric"] in seen:
            return f"Duplicate metric name: {metric['metric']}"
        seen.add(metric["metric"])
    return None


def main() -> int:
    args = parse_args()
    overrides = load_overrides(args.criteria)

    records, invalid, seen_methods, ungated = [], [], set(), []
    for path in sorted(args.results.rglob("*.json")):
        if path.resolve() == args.summary.resolve():
            continue
        try:
            record = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as error:
            invalid.append({"path": str(path), "error": str(error)})
            continue
        error = validate(record)
        if error is None and record["method"] in seen_methods:
            error = f"Duplicate method record: {record['method']}"
        if error:
            invalid.append({"path": str(path), "error": error})
            continue

        for metric in record["metrics"]:
            # criteria.py decides what gates. Whatever the record says about thresholds is
            # a measurement-time artefact and is overwritten here.
            resolved = criterion_for(record["method"], metric["metric"])
            override = overrides.get(f"{record['method']}.{metric['metric']}")
            if override:
                comparison = override.get("comparison") or (resolved[0] if resolved else ">=")
                resolved = (comparison, float(override["tolerance"]), override.get("basis", ""))
                metric["overridden"] = True
            if resolved is None:
                metric.update(
                    comparison="observed", tolerance=None, criterion="recorded, not gating", gating=False, passed=True
                )
                ungated.append(f"{record['method']}.{metric['metric']}")
                continue
            comparison, tolerance, basis = resolved
            metric.update(
                comparison=comparison,
                tolerance=float(tolerance),
                criterion=f"{comparison} {tolerance}",
                gating=True,
                basis=basis,
            )
            metric["passed"] = decide(metric)
            if not metric["passed"]:
                metric["diagnosis"] = diagnosis_for(record["method"], metric["metric"])
        if not any(m["gating"] for m in record["metrics"]):
            invalid.append({"path": record["method"], "error": "No metric in this record has a criterion"})
            continue
        record["passed"] = all(m["passed"] for m in record["metrics"] if m["gating"])

        seen_methods.add(record["method"])
        records.append(record)

    metrics = [m for r in records for m in r["metrics"]]
    gating = [m for m in metrics if bool(m.get("gating", True))]
    summary = {
        "passed": bool(records) and not invalid and all(r["passed"] for r in records),
        "n_methods": len(records),
        "n_passed_methods": sum(bool(r["passed"]) for r in records),
        "n_metrics": len(gating),
        "n_passed_metrics": sum(bool(m["passed"]) for m in gating),
        "n_recorded_metrics": len(metrics),
        "n_informational_metrics": len(metrics) - len(gating),
        "criteria_overrides": sorted(overrides),
        "recorded_without_criterion": ungated,
        "invalid_records": invalid,
        "records": records,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"Evaluated {len(records)} method records, {summary['n_passed_metrics']}/{summary['n_metrics']}"
        f" gating metrics passing -> {args.summary}"
    )
    if overrides:
        print(f"### criteria overridden for: {', '.join(sorted(overrides))}")
    if ungated:
        # Recorded as evidence by design, but listed so a new metric cannot slip in ungated.
        print(f"### {len(ungated)} measurement(s) recorded without a criterion: {', '.join(ungated)}")

    if args.report_dir and args.execution:
        subprocess.run(
            [
                sys.executable,
                str(HERE / "render_results.py"),
                "--summary",
                str(args.summary),
                "--execution",
                str(args.execution),
                "--output-dir",
                str(args.report_dir),
            ],
            check=False,
        )

    if not records:
        raise SystemExit("No result records found. Run the comparison scripts first.")
    if invalid:
        for item in invalid:
            print(f"### invalid record {item['path']}: {item['error']}", file=sys.stderr)
        raise SystemExit(f"Could not parse {len(invalid)} result record(s).")
    if not summary["passed"]:
        raise SystemExit("One or more equivalence criteria failed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
