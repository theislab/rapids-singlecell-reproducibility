"""Fast GPU smoke check to run before the equivalence suite.

Two observed failure modes are infrastructure rather than science, and both waste a
whole run when they are only discovered part-way through the suite. The pinned
`rapids-singlecell-cu12` wheel does not carry every architecture: on V100 (sm_70) CUDA
initializes and then preprocessing kernels fail with `named symbol not found`. A GPU can
also be listed by `nvidia-smi` while CUDA reports `cudaErrorNoDevice`.

This check initializes CuPy, records the hardware and toolkit it actually sees, and
runs one small rapids-singlecell kernel. It always writes its findings as JSON so a
failure leaves durable evidence, and exits non-zero when the GPU cannot run the
suite.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import traceback
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

EXIT_UNUSABLE_GPU = 90


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="Path for the JSON smoke-check record.")
    return parser.parse_args()


# CuPy and rapids-singlecell are installed as CUDA-suffixed wheels, so the plain
# distribution name is not the one that is present.
ALIASES = {
    "rapids-singlecell": ("rapids-singlecell", "rapids-singlecell-cu12"),
    "cupy": ("cupy", "cupy-cuda12x", "cupy-cuda11x"),
}


def package_version(package: str) -> str:
    for candidate in ALIASES.get(package, (package,)):
        try:
            return version(candidate)
        except PackageNotFoundError:
            pass
    return "unknown"


def collect_device_facts() -> dict:
    """Record GPU, driver, CUDA and compute-capability details from CuPy."""
    import cupy as cp

    runtime = cp.cuda.runtime
    device = cp.cuda.Device()
    major, minor = device.compute_capability[0], device.compute_capability[1:]
    attributes = runtime.getDeviceProperties(device.id)
    name = attributes["name"]
    return {
        "device_count": runtime.getDeviceCount(),
        "device_id": device.id,
        "device_name": name.decode() if isinstance(name, bytes) else str(name),
        "compute_capability": f"{major}.{minor}",
        "total_memory_bytes": int(device.mem_info[1]),
        "driver_version": runtime.driverGetVersion(),
        "cuda_runtime_version": runtime.runtimeGetVersion(),
        "cupy_version": package_version("cupy"),
        "cupy_cuda_build": getattr(cp.cuda, "get_local_runtime_version", lambda: None)(),
    }


def run_kernel() -> dict:
    """Exercise one real rapids-singlecell kernel on a tiny matrix.

    `normalize_total` is deliberate: it is the first preprocessing kernel the suite
    touches, and it is where the V100 `named symbol not found` failure surfaced.
    """
    import numpy as np
    import rapids_singlecell as rsc
    from anndata import AnnData

    rng = np.random.default_rng(0)
    counts = rng.poisson(5.0, size=(64, 32)).astype(np.float32)
    adata = AnnData(counts)
    rsc.get.anndata_to_GPU(adata)
    rsc.pp.normalize_total(adata, target_sum=10_000)
    rsc.get.anndata_to_CPU(adata)
    row_sums = np.asarray(adata.X).sum(axis=1)
    return {
        "kernel": "rapids_singlecell.pp.normalize_total",
        "shape": list(adata.shape),
        "max_row_sum_error": float(np.abs(row_sums - 10_000).max()),
    }


def main() -> int:
    args = parse_args()
    record = {
        "checked_at": datetime.now(UTC).isoformat(),
        "hostname": platform.node(),
        "python": platform.python_version(),
        "versions": {
            "rapids-singlecell": package_version("rapids-singlecell"),
            "cupy": package_version("cupy"),
        },
        "passed": False,
        "stage": "import",
        "failure": None,
        "device": None,
        "kernel": None,
    }

    try:
        record["stage"] = "cupy_init"
        record["device"] = collect_device_facts()
        record["stage"] = "rapids_singlecell_kernel"
        record["kernel"] = run_kernel()
        record["passed"] = True
        record["stage"] = "complete"
    except BaseException as error:  # noqa: BLE001 - every failure mode must be recorded
        record["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
            "traceback": traceback.format_exc(),
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2, default=str) + "\n")
    print(json.dumps(record, indent=2, default=str))

    if record["passed"]:
        device = record["device"]
        print(
            f"### gpu smoke check ok: {device['device_name']} "
            f"sm_{device['compute_capability'].replace('.', '')} "
            f"driver={device['driver_version']} cuda={device['cuda_runtime_version']}"
        )
        return 0

    print(f"### gpu smoke check FAILED at stage={record['stage']}", file=sys.stderr)
    print(f"### details persisted to {args.output}", file=sys.stderr)
    return EXIT_UNUSABLE_GPU


if __name__ == "__main__":
    raise SystemExit(main())
