"""Score stored measurements against criteria.py; needs no GPU."""

from __future__ import annotations

import fnmatch
import json
import math
import operator

from benchmarks.comparison import render_results
from benchmarks.comparison.arrays import OUT, captured_points, load_arrays
from benchmarks.comparison.comparisons import AGGREGATIONS, COMPARISONS
from benchmarks.comparison.criteria import CRITERIA, EVIDENCE, criterion_for

COMPARATORS = {"<=": operator.le, "<": operator.lt, ">=": operator.ge, ">": operator.gt}


def decide(metric: dict) -> bool:
    """Apply the metric's criterion; a non-finite measurement fails."""
    if not math.isfinite(float(metric["observed"])):
        return False
    return COMPARATORS[metric["comparison"]](float(metric["observed"]), float(metric["tolerance"]))


def derive(method: str) -> dict[str, float]:
    """Every comparison criteria.py asks for that the stored arrays of one method group can answer."""
    patterns = [rule[0] for rule in CRITERIA.get(method, [])] + EVIDENCE.get(method, [])
    # `markers.mean_set_jaccard` reduces `markers.<group>.set_jaccard`, so those are computed too.
    aggregations = {}
    for pattern in patterns:
        head, _, suffix = pattern.rpartition(".")
        if suffix in AGGREGATIONS:
            member, reduce = AGGREGATIONS[suffix]
            aggregations[pattern] = (f"{head}.*.{member}", reduce)
    wanted = patterns + [member for member, _ in aggregations.values()]

    derived = {}
    for point in captured_points(method):
        for suffix, comparison in COMPARISONS.items():
            name = f"{point}.{suffix}"
            if any(fnmatch.fnmatchcase(name, pattern) for pattern in wanted):
                arrays = load_arrays(method, point, comparison.inputs)
                if arrays is not None:
                    derived[name] = comparison.fn(*arrays)
    for pattern, (member, reduce) in aggregations.items():
        values = [value for name, value in derived.items() if fnmatch.fnmatchcase(name, member)]
        if values:
            derived[pattern] = float(reduce(values))
    return derived


def score(record: dict) -> dict:
    record["metrics"] += [{"metric": name, "observed": value} for name, value in derive(record["method"]).items()]
    # A criterion that matched no measurement fails rather than disappearing.
    record["metrics"] += [
        {"metric": pattern, "observed": math.nan}
        for pattern, *_ in CRITERIA.get(record["method"], [])
        if not any(fnmatch.fnmatchcase(m["metric"], pattern) for m in record["metrics"])
    ]
    for metric in record["metrics"]:
        rule = criterion_for(record["method"], metric["metric"])
        if rule is None:
            metric.update(gating=False, passed=True, comparison="observed", tolerance=None, criterion="recorded")
            continue
        comparison, tolerance = rule
        metric.update(gating=True, comparison=comparison, tolerance=tolerance, criterion=f"{comparison} {tolerance}")
        metric["passed"] = decide(metric)
    gating = [m for m in record["metrics"] if m["gating"]]
    if not gating:
        raise SystemExit(f"No criterion in criteria.py applies to {record['method']}")
    record["passed"] = all(m["passed"] for m in gating)
    return record


def main() -> None:
    records = [score(json.loads(path.read_text())) for path in sorted((OUT / "results").glob("*.json"))]
    if not records:
        raise SystemExit(f"No result records in {OUT / 'results'}. Run the comparison scripts first.")
    gating = [m for r in records for m in r["metrics"] if m["gating"]]
    summary = {
        "passed": all(r["passed"] for r in records),
        "n_metrics": len(gating),
        "n_passed_metrics": sum(m["passed"] for m in gating),
        "records": records,
    }
    (OUT / "equivalence.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Evaluated {len(records)} method records, {summary['n_passed_metrics']}/{summary['n_metrics']} passing")

    # A measurement that is neither gated nor declared evidence should be one or the other.
    undeclared = [
        f"{r['method']}.{m['metric']}"
        for r in records
        for m in r["metrics"]
        if not m["gating"] and not any(fnmatch.fnmatchcase(m["metric"], p) for p in EVIDENCE.get(r["method"], []))
    ]
    for name in undeclared:
        print(f"### recorded without a criterion: {name}")

    execution = OUT / "execution.json"
    if execution.exists():
        render_results.render(summary, json.loads(execution.read_text()), OUT / "report", OUT / "hardware.json")
    if not summary["passed"]:
        raise SystemExit("One or more equivalence criteria failed.")


if __name__ == "__main__":
    main()
