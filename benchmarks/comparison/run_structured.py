from __future__ import annotations

import json
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

from benchmarks.comparison import evaluate
from benchmarks.comparison.arrays import OUT

HERE = Path(__file__).parent
RESULTS = OUT / "results"
LOGS = OUT / "report" / "logs"
SCRIPTS = [
    "squidpy/spatial_autocorr.py",
    "squidpy/co_occurrence.py",
    "squidpy/ligrec.py",
    "scanpy_core/preprocessing.py",
    "scanpy_core/hvg_pca.py",
    "scanpy_core/graphs_embeddings.py",
    "scanpy_core/harmony.py",
    "scanpy_extended/bbknn_scrublet.py",
    "scanpy_extended/embeddings.py",
    "scanpy_extended/ingest_cell_cycle.py",
    "scanpy_extended/rank_genes_groups.py",
    "decoupler/activity.py",
    "pertpy/distance.py",
    "pertpy/guide_assignment.py",
    "pertpy/mixscape.py",
    "pertpy/mixscale.py",
    "biological_pipeline/pbmc3k.py",
]


def record_hardware(output: Path) -> None:
    """Record the GPU and CUDA versions of this run."""
    output.unlink(missing_ok=True)
    try:
        import cupy as cp

        runtime = cp.cuda.runtime
        device = cp.cuda.Device()
        properties = runtime.getDeviceProperties(device.id)
        name = properties["name"]
        facts = {
            "device_name": name.decode() if isinstance(name, bytes) else str(name),
            "compute_capability": f"{properties['major']}.{properties['minor']}",
            "total_memory_bytes": int(device.mem_info[1]),
            "driver_version": runtime.driverGetVersion(),
            "cuda_runtime_version": runtime.runtimeGetVersion(),
        }
    except Exception as error:  # noqa: BLE001 - provenance must never cost the run
        print(f"### hardware not recorded: {type(error).__name__}: {error}", flush=True)
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"measured_at": datetime.now(UTC).isoformat(), "device": facts}, indent=2) + "\n")


requested = sys.argv[1:]
unknown = sorted(set(requested) - set(SCRIPTS))
if unknown:
    raise SystemExit(f"Unknown comparison script(s): {', '.join(unknown)}")
selected = [script for script in SCRIPTS if script in requested] or SCRIPTS
partial = selected != SCRIPTS

RESULTS.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)
if not partial:
    # A partial run keeps the other records, so it still aggregates against the whole suite.
    for previous in [*RESULTS.glob("*.json"), *LOGS.glob("*.log")]:
        previous.unlink()
record_hardware(OUT / "hardware.json")
executions = []

for relative_path in selected:
    script = HERE / relative_path
    log_path = LOGS / f"{relative_path.replace('/', '__')}.log"
    print(f"\n=== {relative_path} ===", flush=True)
    started_at = datetime.now(UTC)
    started = time.monotonic()
    wall = time.time()
    with log_path.open("w") as log_stream:
        completed = subprocess.run(
            [sys.executable, str(script)],
            stdout=log_stream,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    duration = time.monotonic() - started
    log_text = log_path.read_text()
    print(log_text, end="" if log_text.endswith("\n") else "\n")
    result_present = any(path.stat().st_mtime >= wall for path in RESULTS.glob("*.json"))
    passed = completed.returncode == 0 and result_present
    execution = {
        "script": relative_path,
        "exit_code": completed.returncode,
        "passed": passed,
        "started_at": started_at.isoformat(),
        "duration_seconds": duration,
        "log": str(log_path.relative_to(OUT / "report")),
        "result_present": result_present,
    }
    executions.append(execution)

(OUT / "execution.json").write_text(json.dumps({"partial": partial, "executions": executions}, indent=2) + "\n")

# A script without a record is an infrastructure problem, distinct from a criterion being missed.
failed = [execution["script"] for execution in executions if not execution["passed"]]
if failed:
    print(f"### did not produce a record: {', '.join(failed)}")
evaluate.main()
if failed:
    raise SystemExit(f"Comparison scripts did not produce a record: {', '.join(failed)}")
print(f"\nAll comparisons passed. Report: {OUT / 'report' / 'summary.md'}")
