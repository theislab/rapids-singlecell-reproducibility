from __future__ import annotations

import json
import os
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any


def package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "unknown"


def upper_bound(metric: str, observed: float, threshold: float) -> dict[str, Any]:
    return {
        "metric": metric,
        "observed": float(observed),
        "comparison": "<=",
        "tolerance": float(threshold),
        "criterion": f"<= {threshold}",
        "passed": bool(observed <= threshold),
    }


def lower_bound(metric: str, observed: float, threshold: float) -> dict[str, Any]:
    return {
        "metric": metric,
        "observed": float(observed),
        "comparison": ">=",
        "tolerance": float(threshold),
        "criterion": f">= {threshold}",
        "passed": bool(observed >= threshold),
    }


def write_report(
    *, method: str, dataset: str, tier: str, metrics: list[dict[str, Any]]
) -> Path:
    report = {
        "method": method,
        "reference_package": "squidpy",
        "dataset": dataset,
        "tier": tier,
        "versions": {
            "rapids-singlecell": package_version("rapids-singlecell"),
            "squidpy": package_version("squidpy"),
        },
        "passed": all(metric["passed"] for metric in metrics),
        "metrics": metrics,
    }

    output_dir = Path(
        os.environ.get(
            "EQUIVALENCE_OUTPUT_DIR", Path(__file__).parent / "results"
        )
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{method}.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))

    failures = [metric for metric in metrics if not metric["passed"]]
    if failures:
        details = ", ".join(
            f"{item['metric']}={item['observed']} ({item['criterion']})"
            for item in failures
        )
        raise AssertionError(f"Equivalence thresholds failed: {details}")
    return output
