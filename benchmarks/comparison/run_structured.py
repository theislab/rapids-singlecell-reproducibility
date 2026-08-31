from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from arrays import OUT

RESULTS = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", OUT / "results"))
SUMMARY = Path(os.environ.get("EQUIVALENCE_SUMMARY", OUT / "equivalence.json"))
EXECUTION = Path(os.environ.get("EQUIVALENCE_EXECUTION", OUT / "execution.json"))
REPORT_DIR = Path(os.environ.get("EQUIVALENCE_REPORT_DIR", OUT / "report"))
ARRAYS = Path(os.environ.get("EQUIVALENCE_ARRAY_DIR", OUT / "arrays"))
HARDWARE = REPORT_DIR.parent / "hardware.json"
SCRIPTS = [
    # The expensive niche comparison runs last.
    "squidpy/spatial_autocorr.py",
    "squidpy/co_occurrence.py",
    "squidpy/ligrec.py",
    "scanpy_core/preprocessing.py",
    "scanpy_core/hvg_pca.py",
    "scanpy_core/graphs_embeddings.py",
    "scanpy_core/harmony.py",
    "scanpy_extended/sqrt.py",
    "scanpy_extended/bbknn_scrublet.py",
    "scanpy_extended/clustering.py",
    "scanpy_extended/embeddings.py",
    "scanpy_extended/ingest_cell_cycle.py",
    "scanpy_extended/rank_genes_groups.py",
    "decoupler/decoupler.py",
    "pertpy/distance.py",
    "pertpy/guide_assignment.py",
    "pertpy/mixscape.py",
    "pertpy/mixscale.py",
    "biological_pipeline/pbmc3k.py",
    "squidpy/calculate_niche.py",
]


def record_hardware(output: Path) -> None:
    """Name the GPU and CUDA versions this run measured on. The report travels without the
    run, so a number in it is only interpretable next to the hardware that produced it."""
    try:
        import cupy as cp

        runtime = cp.cuda.runtime
        device = cp.cuda.Device()
        major, minor = device.compute_capability[0], device.compute_capability[1:]
        name = runtime.getDeviceProperties(device.id)["name"]
        facts = {
            "device_name": name.decode() if isinstance(name, bytes) else str(name),
            "compute_capability": f"{major}.{minor}",
            "total_memory_bytes": int(device.mem_info[1]),
            "driver_version": runtime.driverGetVersion(),
            "cuda_runtime_version": runtime.runtimeGetVersion(),
        }
    except Exception as error:  # noqa: BLE001 - provenance must never cost the run
        print(f"### hardware not recorded: {type(error).__name__}: {error}", flush=True)
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"measured_at": datetime.now(UTC).isoformat(), "device": facts}, indent=2) + "\n")


def select_scripts(argv: list[str]) -> tuple[list[str], bool]:
    """Resolve which comparisons to run, from arguments or EQUIVALENCE_SCRIPTS.

    A subset run is for iterating on one method; it is reported as partial so its
    aggregate output is never mistaken for a full-suite run.
    """
    requested = list(argv) or [
        item for item in os.environ.get("EQUIVALENCE_SCRIPTS", "").replace(",", " ").split() if item
    ]
    if not requested:
        return SCRIPTS, False
    unknown = [item for item in requested if item not in SCRIPTS]
    if unknown:
        raise SystemExit(f"Unknown comparison script(s): {', '.join(unknown)}\nAvailable:\n  " + "\n  ".join(SCRIPTS))
    # Keep the curated order regardless of the order given on the command line.
    return [item for item in SCRIPTS if item in set(requested)], True


selected, partial = select_scripts(sys.argv[1:])

RESULTS.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)
log_dir = REPORT_DIR / "logs"
log_dir.mkdir(parents=True, exist_ok=True)
if partial:
    print(f"### partial run: {len(selected)}/{len(SCRIPTS)} comparisons", flush=True)
    # Records from comparisons that are not rerun are left in place on purpose, so a
    # targeted rerun still aggregates against the rest of the suite.
    for relative_path in selected:
        stale_log = log_dir / f"{relative_path.replace('/', '__')}.log"
        stale_log.unlink(missing_ok=True)
else:
    for previous_result in RESULTS.glob("*.json"):
        previous_result.unlink()
    for previous_log in log_dir.glob("*.log"):
        previous_log.unlink()

record_hardware(HARDWARE)

env = os.environ.copy()
env["EQUIVALENCE_OUTPUT_DIR"] = str(RESULTS)
env["EQUIVALENCE_ARRAY_DIR"] = str(ARRAYS)
env.setdefault("EQUIVALENCE_ARTIFACT_DIR", str(REPORT_DIR / "artifacts"))
failures = []
executions = []

for relative_path in selected:
    script = HERE / relative_path
    log_path = log_dir / f"{relative_path.replace('/', '__')}.log"
    print(f"\n=== {relative_path} ===", flush=True)
    started_at = datetime.now(UTC)
    started = time.monotonic()
    previous_results = {path.resolve(): path.stat().st_mtime_ns for path in RESULTS.glob("*.json")}
    with log_path.open("w") as log_stream:
        completed = subprocess.run(
            [sys.executable, str(script)],
            env=env,
            stdout=log_stream,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    duration = time.monotonic() - started
    log_text = log_path.read_text()
    print(log_text, end="" if log_text.endswith("\n") else "\n")
    result_present = any(
        path.resolve() not in previous_results or path.stat().st_mtime_ns != previous_results[path.resolve()]
        for path in RESULTS.glob("*.json")
    )
    passed = completed.returncode == 0 and result_present
    execution = {
        "script": relative_path,
        "exit_code": completed.returncode,
        "passed": passed,
        "started_at": started_at.isoformat(),
        "duration_seconds": duration,
        "log": str(log_path.relative_to(REPORT_DIR)),
        "result_present": result_present,
    }
    executions.append(execution)
    if not passed:
        failures.append(
            f"{relative_path} (exit {completed.returncode}, result record {'present' if result_present else 'missing'})"
        )

EXECUTION.parent.mkdir(parents=True, exist_ok=True)
EXECUTION.write_text(
    json.dumps(
        {
            "partial": partial,
            "n_selected": len(selected),
            "n_available": len(SCRIPTS),
            "executions": executions,
        },
        indent=2,
    )
    + "\n"
)

# `failures` above are scripts that produced no record at all: an infrastructure problem,
# distinct from a criterion being missed, which evaluate.py decides.
evaluation = subprocess.run(
    [
        sys.executable,
        str(HERE / "evaluate.py"),
        "--results",
        str(RESULTS),
        "--summary",
        str(SUMMARY),
        "--execution",
        str(EXECUTION),
        "--report-dir",
        str(REPORT_DIR),
        "--arrays",
        str(ARRAYS),
        "--write-enriched",
    ],
    env=env,
    check=False,
)

if failures:
    raise SystemExit(f"Comparison scripts did not produce a record: {', '.join(failures)}")
if evaluation.returncode:
    raise SystemExit(f"Evaluation failed with exit {evaluation.returncode}")

print(f"\nAll structured comparisons passed. Summary: {SUMMARY}; report: {REPORT_DIR / 'summary.md'}")
