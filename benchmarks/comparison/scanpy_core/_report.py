from __future__ import annotations

import json
import os
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def package_version(package: str) -> str:
    candidates = (package, "rapids-singlecell-cu12") if package == "rapids-singlecell" else (package,)
    for candidate in candidates:
        try:
            return version(candidate)
        except PackageNotFoundError:
            pass
    return "unknown"


def upper(metric: str, observed: float, threshold: float) -> dict:
    return {
        "metric": metric,
        "observed": float(observed),
        "comparison": "<=",
        "tolerance": float(threshold),
        "criterion": f"<= {threshold}",
        "passed": bool(observed <= threshold),
    }


def lower(metric: str, observed: float, threshold: float) -> dict:
    return {
        "metric": metric,
        "observed": float(observed),
        "comparison": ">=",
        "tolerance": float(threshold),
        "criterion": f">= {threshold}",
        "passed": bool(observed >= threshold),
    }


def write_report(
    method: str,
    dataset: str,
    tier: str,
    metrics: list[dict],
    *,
    reference_package: str = "scanpy",
    packages: tuple[str, ...] = ("scanpy",),
) -> Path:
    report = {
        "method": method,
        "reference_package": reference_package,
        "dataset": dataset,
        "tier": tier,
        "versions": {
            "rapids-singlecell": package_version("rapids-singlecell"),
            **{package: package_version(package) for package in packages},
        },
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
        raise AssertionError(
            "Equivalence thresholds failed: "
            + ", ".join(f"{metric['metric']}={metric['observed']} ({metric['criterion']})" for metric in failed)
        )
    return output
