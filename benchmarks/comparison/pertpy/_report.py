from __future__ import annotations

import json
import os
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def package_version(package):
    candidates = (package, "rapids-singlecell-cu12") if package == "rapids-singlecell" else (package,)
    for candidate in candidates:
        try:
            return version(candidate)
        except PackageNotFoundError:
            pass
    return "unknown"


def upper(name, value, threshold):
    return {
        "metric": name,
        "observed": float(value),
        "comparison": "<=",
        "tolerance": float(threshold),
        "criterion": f"<= {threshold}",
        "passed": bool(value <= threshold),
    }


def lower(name, value, threshold):
    return {
        "metric": name,
        "observed": float(value),
        "comparison": ">=",
        "tolerance": float(threshold),
        "criterion": f">= {threshold}",
        "passed": bool(value >= threshold),
    }


def write_report(method, dataset, tier, metrics):
    report = {
        "method": method,
        "reference_package": "pertpy",
        "dataset": dataset,
        "tier": tier,
        "versions": {"rapids-singlecell": package_version("rapids-singlecell"), "pertpy": package_version("pertpy")},
        "passed": all(item["passed"] for item in metrics),
        "metrics": metrics,
    }
    output_dir = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", Path(__file__).parent / "results"))
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / f"{method}.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise AssertionError(f"{method} equivalence thresholds failed")
