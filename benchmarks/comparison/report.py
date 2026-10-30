"""Write one method group's measurement record."""

from __future__ import annotations

import json
from importlib.metadata import version

from benchmarks.comparison.arrays import OUT


def measure(metric: str, observed: float) -> dict:
    return {"metric": metric, "observed": float(observed)}


def write_report(
    method: str,
    dataset: str,
    metrics: list[dict],
    *,
    reference_package: str = "scanpy",
    shape: tuple[int, int] | None = None,
) -> None:
    """Store one group's record; `shape` is the compared input as `(cells, genes)`."""
    report = {
        "method": method,
        "reference_package": reference_package,
        "dataset": dataset,
        "versions": {
            "rapids-singlecell": version("rapids-singlecell-cu13"),
            reference_package: version(reference_package),
        },
        "metrics": metrics,
        **({"shape": list(shape)} if shape is not None else {}),
    }
    output = OUT / "results" / f"{method}.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
