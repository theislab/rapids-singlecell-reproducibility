from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).parent
RESULTS = HERE / "results"
SCRIPTS = [
    "scanpy_extended/sqrt.py",
    "scanpy_extended/bbknn_scrublet.py",
    "scanpy_extended/clustering.py",
    "scanpy_extended/embeddings.py",
    "scanpy_extended/ingest_cell_cycle.py",
    "scanpy_extended/rank_genes_groups.py",
    "squidpy/spatial_autocorr.py",
    "squidpy/co_occurrence.py",
    "squidpy/ligrec.py",
    "squidpy/calculate_niche.py",
    "decoupler/decoupler.py",
    "pertpy/distance.py",
    "pertpy/guide_assignment.py",
    "pertpy/mixscape.py",
    "pertpy/mixscale.py",
]

RESULTS.mkdir(parents=True, exist_ok=True)
for previous_result in RESULTS.glob("*.json"):
    previous_result.unlink()
env = os.environ.copy()
env["EQUIVALENCE_OUTPUT_DIR"] = str(RESULTS)
failures = []

for relative_path in SCRIPTS:
    script = HERE / relative_path
    print(f"\n=== {relative_path} ===", flush=True)
    completed = subprocess.run([sys.executable, str(script)], env=env, check=False)
    if completed.returncode:
        failures.append((relative_path, completed.returncode))

summary = HERE / "equivalence.json"
env["EQUIVALENCE_SUMMARY"] = str(summary)
aggregation = subprocess.run(
    [sys.executable, str(HERE / "collect_results.py")], env=env, check=False
)

if failures:
    detail = ", ".join(f"{path} (exit {code})" for path, code in failures)
    raise SystemExit(f"Comparison failures: {detail}")
if aggregation.returncode:
    raise SystemExit(f"Result aggregation failed with exit {aggregation.returncode}")

print(f"\nAll structured comparisons passed. Summary: {summary}")
