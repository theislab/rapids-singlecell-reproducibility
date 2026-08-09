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


def measure(metric: str, observed: float) -> dict:
    """Record one measurement. Nothing here decides anything.

    Whether this metric gates, and against what threshold, is `criteria.py`'s business and
    `evaluate.py` applies it. Keeping the two apart is what lets a stored run be re-scored
    without a GPU.
    """
    return {"metric": metric, "observed": float(observed)}


def write_report(method, dataset, tier, metrics):
    report = {
        "method": method,
        "reference_package": "pertpy",
        "dataset": dataset,
        "tier": tier,
        "versions": {"rapids-singlecell": package_version("rapids-singlecell"), "pertpy": package_version("pertpy")},
        "metrics": metrics,
    }
    output_dir = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", Path(__file__).parent / "results"))
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / f"{method}.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    # Measuring is not evaluating: the verdict is formed by evaluate.py from this
    # record, so a missed threshold here is evidence, not a script failure.
