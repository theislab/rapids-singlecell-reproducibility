"""Write one method group's measurement record. Nothing here decides anything.

Whether a metric gates, and against what threshold, is `criteria.py`'s business and
`evaluate.py` applies it. Keeping the two apart is what lets a stored run be re-scored
without a GPU, so a missed threshold at measurement time is evidence, not a script failure.

Each group directory keeps a `_report.py` that re-exports this module and binds the
reference package that group compares against.
"""

from __future__ import annotations

import json
import os
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from arrays import OUT, capture  # capture is re-exported; the comparison scripts import it from here


def package_version(package: str) -> str:
    candidates = (package, "rapids-singlecell-cu12") if package == "rapids-singlecell" else (package,)
    for candidate in candidates:
        try:
            return version(candidate)
        except PackageNotFoundError:
            pass
    return "unknown"


def measure(metric: str, observed: float) -> dict:
    return {"metric": metric, "observed": float(observed)}


def write_report(
    method: str,
    dataset: str,
    tier: str,
    metrics: list[dict],
    *,
    reference_package: str = "scanpy",
    packages: tuple[str, ...] | None = None,
    shape: tuple[int, int] | None = None,
) -> Path:
    """Store one group's record. `shape` is the compared input as `(cells, genes)`, so the
    report states the scale each comparison ran at instead of it being written by hand."""
    report = {
        "method": method,
        "reference_package": reference_package,
        "dataset": dataset,
        "tier": tier,
        "versions": {
            "rapids-singlecell": package_version("rapids-singlecell"),
            **{package: package_version(package) for package in packages or (reference_package,)},
        },
        "metrics": metrics,
        **({"shape": list(shape)} if shape is not None else {}),
    }
    output_dir = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", OUT / "results"))
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{method}.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return output
