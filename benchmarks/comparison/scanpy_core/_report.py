from __future__ import annotations

import json
import os
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1]))

from arrays import capture


def package_version(package: str) -> str:
    candidates = (package, "rapids-singlecell-cu12") if package == "rapids-singlecell" else (package,)
    for candidate in candidates:
        try:
            return version(candidate)
        except PackageNotFoundError:
            pass
    return "unknown"


def measure(metric: str, observed: float) -> dict:
    """Record one measurement. Nothing here decides anything.

    Whether this metric gates, and against what threshold, is `criteria.py`'s business and
    `evaluate.py` applies it. Keeping the two apart is what lets a stored run be re-scored
    without a GPU.
    """
    return {"metric": metric, "observed": float(observed)}


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
        "metrics": metrics,
    }
    output_dir = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", Path(__file__).parent / "results"))
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{method}.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return output
