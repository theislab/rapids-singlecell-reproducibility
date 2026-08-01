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


def upper(metric: str, observed: float, threshold: float, *, basis: str = "") -> dict:
    return {
        "metric": metric,
        "observed": float(observed),
        "comparison": "<=",
        "tolerance": float(threshold),
        "criterion": f"<= {threshold}",
        "passed": bool(observed <= threshold),
        "gating": True,
        "basis": basis,
    }


def lower(metric: str, observed: float, threshold: float, *, basis: str = "") -> dict:
    return {
        "metric": metric,
        "observed": float(observed),
        "comparison": ">=",
        "tolerance": float(threshold),
        "criterion": f">= {threshold}",
        "passed": bool(observed >= threshold),
        "gating": True,
        "basis": basis,
    }


def observed_only(metric: str, observed: float, *, basis: str) -> dict:
    """Record a measurement that is evidence rather than a pass/fail criterion.

    Absolute quality scores of a single implementation (for example UMAP
    trustworthiness) do not test CPU/GPU equivalence, so they are reported
    without gating the suite. The paired equivalence criterion derived from them
    is what gates.
    """
    return {
        "metric": metric,
        "observed": float(observed),
        "comparison": "observed",
        "tolerance": None,
        "criterion": "recorded, not gating",
        "passed": True,
        "gating": False,
        "basis": basis,
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
    gating = [metric for metric in metrics if metric.get("gating", True)]
    report = {
        "method": method,
        "reference_package": reference_package,
        "dataset": dataset,
        "tier": tier,
        "versions": {
            "rapids-singlecell": package_version("rapids-singlecell"),
            **{package: package_version(package) for package in packages},
        },
        "passed": all(metric["passed"] for metric in gating),
        "metrics": metrics,
    }
    output_dir = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", Path(__file__).parent / "results"))
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{method}.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    failed = [metric for metric in gating if not metric["passed"]]
    if failed:
        raise AssertionError(
            "Equivalence thresholds failed: "
            + ", ".join(f"{metric['metric']}={metric['observed']} ({metric['criterion']})" for metric in failed)
        )
    return output
