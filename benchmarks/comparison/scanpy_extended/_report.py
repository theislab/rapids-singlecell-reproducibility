from __future__ import annotations

import json
import os
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def _version(package):
    try:
        return version(package)
    except PackageNotFoundError:
        return "unknown"


def upper(metric, observed, threshold):
    return {"metric": metric, "observed": float(observed), "comparison": "<=", "tolerance": float(threshold), "criterion": f"<= {threshold}", "passed": bool(observed <= threshold)}


def lower(metric, observed, threshold):
    return {"metric": metric, "observed": float(observed), "comparison": ">=", "tolerance": float(threshold), "criterion": f">= {threshold}", "passed": bool(observed >= threshold)}


def write_report(method, dataset, tier, metrics):
    report = {
        "method": method,
        "reference_package": "scanpy",
        "dataset": dataset,
        "tier": tier,
        "versions": {"rapids-singlecell": _version("rapids-singlecell"), "scanpy": _version("scanpy")},
        "passed": all(metric["passed"] for metric in metrics),
        "metrics": metrics,
    }
    output_dir = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", Path(__file__).parent / "results"))
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{method}.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    failed = [metric for metric in metrics if not metric["passed"]]
    if failed:
        raise AssertionError("Equivalence thresholds failed: " + ", ".join(f"{m['metric']}={m['observed']} ({m['criterion']})" for m in failed))
