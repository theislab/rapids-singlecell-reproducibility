"""Evaluate stored measurements against their criteria. No GPU, no rerun. See README.md.

Criteria come from `criteria.py`, never from the record, so an old result directory can be
re-scored under today's criteria.

    python benchmarks/comparison/evaluate.py
    python benchmarks/comparison/evaluate.py --results <previous-run>/results
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
import fnmatch
import json
import math
import operator
import os
import subprocess
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from arrays import OUT, captured_points, load_arrays
from comparisons import AGGREGATIONS, ALLCLOSE_EVIDENCE, COMPARISONS
from criteria import CRITERIA, EVIDENCE, EVIDENCE_BASIS, criterion_for, diagnosis_for

HERE = Path(__file__).parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--results",
        type=Path,
        default=Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", OUT / "results")),
        help="Directory of per-method measurement records.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path(os.environ.get("EQUIVALENCE_SUMMARY", OUT / "equivalence.json")),
        help="Where to write the aggregated evaluation.",
    )
    parser.add_argument("--criteria", type=Path, help='TOML file overriding thresholds, keyed "<method>.<metric>".')
    parser.add_argument(
        "--arrays",
        type=Path,
        help="Directory of per-method Zarr stores holding the raw outputs. Comparisons are computed "
        "from these when present; otherwise the scalar the script recorded is used.",
    )
    parser.add_argument("--execution", type=Path, help="Execution record; enables report rendering.")
    parser.add_argument("--report-dir", type=Path, help="Render the reviewer report here.")
    parser.add_argument(
        "--hardware",
        type=Path,
        help="Hardware record, for the provenance block in the report. "
        "Defaults to hardware.json beside the report directory, where a run writes it.",
    )
    parser.add_argument(
        "--write-enriched",
        action="store_true",
        help="Rewrite each result record with the comparisons derived from the arrays inlined. "
        "Records in this form re-score without the Zarr stores, which are not kept.",
    )
    return parser.parse_args()


def load_overrides(path: Path | None) -> dict:
    if path is None:
        return {}
    overrides = tomllib.loads(path.read_text())
    for metric, entry in overrides.items():
        if "tolerance" not in entry:
            raise SystemExit(f"Criteria override for {metric!r} must set a tolerance")
    return overrides


COMPARATORS = {"<=": operator.le, "<": operator.lt, ">=": operator.ge, ">": operator.gt}


def decide(metric: dict) -> bool:
    """Apply the metric's criterion. A non-finite measurement fails; an unknown
    comparator is a hard error rather than a silently substituted default."""
    if not math.isfinite(float(metric["observed"])):
        return False
    compare = COMPARATORS.get(metric["comparison"])
    if compare is None:
        raise SystemExit(f"Unknown comparison {metric['comparison']!r} for metric {metric['metric']!r}")
    return compare(float(metric["observed"]), float(metric["tolerance"]))


def enriched(record: dict) -> dict:
    """The record with derived comparisons inlined, so it re-scores without the Zarr
    stores. Measurements and the verdict only -- never a copy of the criterion."""
    return {
        **{key: value for key, value in record.items() if key not in ("metrics", "derived_metrics", "passed")},
        "metrics": [{"metric": m["metric"], "observed": m["observed"]} for m in record["metrics"]],
        "passed": record["passed"],
    }


def validate(record: dict) -> str | None:
    """Reject a malformed measurement record rather than silently scoring it.

    Only measurements are required here. Anything a script may still write about
    thresholds is ignored: criteria.py is the authority on what a criterion is.
    """
    if not {"method", "reference_package", "dataset", "tier", "versions", "metrics"}.issubset(record):
        return "Missing required result fields"
    # An empty list is normal now: a script whose comparisons all come from stored arrays
    # records no scalars of its own.
    if not isinstance(record["metrics"], list):
        return "Metrics must be a list"
    seen = set()
    for metric in record["metrics"]:
        if not {"metric", "observed"}.issubset(metric):
            return "Metric is missing required fields"
        if metric["metric"] in seen:
            return f"Duplicate metric name: {metric['metric']}"
        seen.add(metric["metric"])
    return None


MISAPPLIED: list[str] = []


def derive_from_arrays(method: str) -> dict[str, float]:
    """Compute every comparison the stored arrays support for one method group.

    A criterion named `<point>.<suffix>` is computed when `<point>` carries the arrays that
    `<suffix>` declares it needs. Evidence quantities for an `allclose` criterion, and
    aggregations over a family of metrics, are derived alongside — so a question about a
    run never has to be anticipated while the run is happening.
    """
    points = captured_points(method)
    if not points:
        return {}

    wanted: dict[str, str] = {}
    # Evidence patterns are derived exactly like criteria; they simply resolve to no rule
    # later, so they are reported without a verdict.
    for pattern in [rule[0] for rule in CRITERIA.get(method, [])] + EVIDENCE.get(method, []):
        head, _, suffix = pattern.rpartition(".")
        # `markers.mean_set_jaccard` is a reduction over `markers.<group>.set_jaccard`, so
        # the members have to be computed even though no criterion names them directly.
        if suffix in AGGREGATIONS:
            member = AGGREGATIONS[suffix][0]
            for point in points:
                if point.startswith(f"{head}."):
                    wanted[f"{point}.{member}"] = member
            continue
        # Some suffixes are themselves dotted, e.g. `umap.cpu.trustworthiness`.
        for candidate_suffix in (suffix, ".".join(pattern.rsplit(".", 2)[-2:])):
            if candidate_suffix not in COMPARISONS:
                continue
            head = pattern[: -(len(candidate_suffix) + 1)]
            for point in points:
                if head in (point, "*") or fnmatch.fnmatchcase(point, head):
                    wanted[f"{point}.{candidate_suffix}"] = candidate_suffix
                    if candidate_suffix.endswith("allclose_excess"):
                        # Mirror the criterion's own prefix: evidence for an `abs_` criterion
                        # has to be measured on magnitudes too, or a sign flip reads as a
                        # catastrophic error beside a criterion that passed.
                        prefix = candidate_suffix[: -len("allclose_excess")]
                        for evidence in ALLCLOSE_EVIDENCE:
                            variant = f"{prefix}{evidence}"
                            name = variant if variant in COMPARISONS else evidence
                            wanted[f"{point}.{name}"] = name
            break

    derived: dict[str, float] = {}
    for name, suffix in wanted.items():
        point = name[: -(len(suffix) + 1)]
        comparison = COMPARISONS[suffix]
        loaded = load_arrays(method, point, comparison.inputs)
        if loaded is None:
            continue
        try:
            derived[name] = comparison.fn(*loaded)
        except Exception as error:  # noqa: BLE001 - a bad rule must not cost the whole run
            # Almost always a criterion pointed at arrays the comparison cannot handle, e.g.
            # a numeric comparison matched against captured gene names by too broad a glob.
            MISAPPLIED.append(f"{method}.{name}: {type(error).__name__}: {error}")

    # Criteria stated over a family — "the weakest marker overlap across cell types" — are
    # aggregations of metrics just computed, not separate measurements.
    for pattern, *_ in CRITERIA.get(method, []):
        head, _, suffix = pattern.rpartition(".")
        if suffix not in AGGREGATIONS:
            continue
        member_suffix, reduce = AGGREGATIONS[suffix]
        values = [v for k, v in derived.items() if k.endswith(f".{member_suffix}") and k.startswith(f"{head}.")]
        if values:
            derived[pattern] = reduce(values)
    return derived


def main() -> int:
    args = parse_args()
    overrides = load_overrides(args.criteria)

    if args.arrays:
        os.environ["EQUIVALENCE_ARRAY_DIR"] = str(args.arrays)

    records, invalid, seen_methods, ungated, sources = [], [], set(), [], {}
    declared_evidence: set[str] = set()
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

        # Comparisons the stored arrays can answer are computed here rather than trusted
        # from the record, which is what makes a new comparison free of a GPU run.
        for name, observed in derive_from_arrays(record["method"]).items():
            existing = next((m for m in record["metrics"] if m["metric"] == name), None)
            if existing is None:
                record["metrics"].append({"metric": name, "observed": observed})
            else:
                existing["observed"] = observed
            record.setdefault("derived_metrics", []).append(name)

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
                name = f"{record['method']}.{metric['metric']}"
                ungated.append(name)
                if any(fnmatch.fnmatchcase(metric["metric"], p) for p in EVIDENCE.get(record["method"], [])):
                    declared_evidence.add(name)
                    metric["basis"] = EVIDENCE_BASIS
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
            invalid.append(
                {
                    "path": record["method"],
                    "error": "No criterion applies: the script captured nothing criteria.py knows about",
                }
            )
            continue
        record["passed"] = all(m["passed"] for m in record["metrics"] if m["gating"])

        seen_methods.add(record["method"])
        records.append(record)
        sources[record["method"]] = path

    if args.write_enriched:
        for record in records:
            sources[record["method"]].write_text(json.dumps(enriched(record), indent=2) + "\n")

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
        "misapplied_comparisons": MISAPPLIED,
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
        # Evidence derived beside an allclose criterion, and comparisons criteria.py declares
        # as evidence outright, are ungated by design and there is a lot of both; name only
        # the rest, so a metric that should gate cannot hide in the noise.
        expected = tuple(f".{name}" for name in ALLCLOSE_EVIDENCE)
        notable = [name for name in ungated if not name.endswith(expected) and name not in declared_evidence]
        print(
            f"### {len(ungated)} measurement(s) recorded without a criterion, {len(notable)} of them not allclose evidence"
        )
        for name in notable:
            print(f"###   {name}")

    if args.report_dir and args.execution:
        # The report is promoted as a single file and travels without the run, so it names the
        # GPU it was measured on. run_structured.py writes this beside the report directory.
        hardware = args.hardware or args.report_dir.parent / "hardware.json"
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
                *(["--hardware", str(hardware)] if hardware.exists() else []),
            ],
            check=False,
        )

    if MISAPPLIED:
        for item in MISAPPLIED:
            print(f"### comparison could not be applied: {item}", file=sys.stderr)
    if not records:
        raise SystemExit("No result records found. Run the comparison scripts first.")
    if MISAPPLIED:
        raise SystemExit(f"{len(MISAPPLIED)} criterion/criteria point at arrays their comparison cannot use.")
    if invalid:
        for item in invalid:
            print(f"### invalid record {item['path']}: {item['error']}", file=sys.stderr)
        raise SystemExit(f"Could not parse {len(invalid)} result record(s).")
    if not summary["passed"]:
        raise SystemExit("One or more equivalence criteria failed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
