from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

HERE = Path(__file__).parent
RESULTS = Path(os.environ.get("EQUIVALENCE_OUTPUT_DIR", HERE / "results"))
SUMMARY = Path(os.environ.get("EQUIVALENCE_SUMMARY", HERE / "equivalence.json"))
EXECUTION = Path(os.environ.get("EQUIVALENCE_EXECUTION", HERE / "execution.json"))
REPORT_DIR = Path(os.environ.get("EQUIVALENCE_REPORT_DIR", HERE / "report"))
SCRIPTS = [
    # Prioritize the Squidpy comparisons requested by the reviewers. The more
    # expensive niche comparison intentionally runs last.
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


RESULTS.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)
log_dir = REPORT_DIR / "logs"
log_dir.mkdir(parents=True, exist_ok=True)
for previous_result in RESULTS.glob("*.json"):
    previous_result.unlink()
for previous_log in log_dir.glob("*.log"):
    previous_log.unlink()

env = os.environ.copy()
env["EQUIVALENCE_OUTPUT_DIR"] = str(RESULTS)
env.setdefault("EQUIVALENCE_ARTIFACT_DIR", str(REPORT_DIR / "artifacts"))
failures = []
executions = []

for relative_path in SCRIPTS:
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
EXECUTION.write_text(json.dumps({"executions": executions}, indent=2) + "\n")

env["EQUIVALENCE_SUMMARY"] = str(SUMMARY)
aggregation = subprocess.run([sys.executable, str(HERE / "collect_results.py")], env=env, check=False)
rendering = subprocess.run(
    [
        sys.executable,
        str(HERE / "render_results.py"),
        "--summary",
        str(SUMMARY),
        "--execution",
        str(EXECUTION),
        "--output-dir",
        str(REPORT_DIR),
    ],
    check=False,
)

if failures:
    raise SystemExit(f"Comparison failures: {', '.join(failures)}")
if aggregation.returncode:
    raise SystemExit(f"Result aggregation failed with exit {aggregation.returncode}")
if rendering.returncode:
    raise SystemExit(f"Result rendering failed with exit {rendering.returncode}")

print(f"\nAll structured comparisons passed. Summary: {SUMMARY}; report: {REPORT_DIR / 'summary.md'}")
