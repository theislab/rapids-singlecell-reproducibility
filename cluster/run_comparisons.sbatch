#!/usr/bin/env bash
#SBATCH --job-name=rsc-equivalence
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_priority
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=90G
#SBATCH --time=04:00:00
#SBATCH --nice=10000
#SBATCH --output=rsc-equivalence-%j.out
#SBATCH --error=rsc-equivalence-%j.err

set -euo pipefail

# The comparison repository lives on shared storage and is copied to node-local
# scratch before any imports, package installs, or comparisons run.
REPRO_SOURCE="${REPRO_SOURCE:-$SLURM_SUBMIT_DIR}"
DURABLE_OUTPUT="${EQUIVALENCE_DURABLE_OUTPUT:-$REPRO_SOURCE/benchmarks/comparison/cluster-results/$SLURM_JOB_ID}"

LOCAL_BASE="${SLURM_TMPDIR:-/localscratch/$USER}/rapids-singlecell-equivalence"
JOB_ROOT="$LOCAL_BASE/jobs/$SLURM_JOB_ID"
CACHE_ROOT="$LOCAL_BASE/cache"
LOCAL_REPRO="$JOB_ROOT/rapids-singlecell-reproducibility"

mkdir -p "$JOB_ROOT" "$CACHE_ROOT" "$DURABLE_OUTPUT"
mkdir -p "$CACHE_ROOT/uv" "$CACHE_ROOT/xdg" "$CACHE_ROOT/cupy" "$CACHE_ROOT/numba"
mkdir -p "$CACHE_ROOT/matplotlib" "$CACHE_ROOT/tmp/$SLURM_JOB_ID" "$CACHE_ROOT/joblib/$SLURM_JOB_ID"

export UV_CACHE_DIR="$CACHE_ROOT/uv"
export UV_PROJECT_ENVIRONMENT="$JOB_ROOT/.venv"
export UV_COMPILE_BYTECODE=1
export UV_PYTHON=3.12
export XDG_CACHE_HOME="$CACHE_ROOT/xdg"
export CUPY_CACHE_DIR="$CACHE_ROOT/cupy"
export NUMBA_CACHE_DIR="$CACHE_ROOT/numba"
export MPLCONFIGDIR="$CACHE_ROOT/matplotlib"
export TMPDIR="$CACHE_ROOT/tmp/$SLURM_JOB_ID"
export JOBLIB_TEMP_FOLDER="$CACHE_ROOT/joblib/$SLURM_JOB_ID"
export EQUIVALENCE_OUTPUT_DIR="$LOCAL_REPRO/benchmarks/comparison/results"
export EQUIVALENCE_SUMMARY="$LOCAL_REPRO/benchmarks/comparison/equivalence.json"
export EQUIVALENCE_EXECUTION="$LOCAL_REPRO/benchmarks/comparison/execution.json"
export EQUIVALENCE_REPORT_DIR="$LOCAL_REPRO/benchmarks/comparison/report"
export EQUIVALENCE_ARTIFACT_DIR="$EQUIVALENCE_REPORT_DIR/artifacts"
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

echo "### host=$(hostname) job=$SLURM_JOB_ID started=$(date --iso-8601=seconds)"
echo "### source=$REPRO_SOURCE"
echo "### local_base=$LOCAL_BASE"
echo "### durable_output=$DURABLE_OUTPUT"
nvidia-smi

rsync -a --delete --exclude=.git --exclude=.venv --exclude=cluster-results \
    "$REPRO_SOURCE/" "$LOCAL_REPRO/"

cd "$LOCAL_REPRO"
if ! command -v uv >/dev/null; then
    UV_TOOL_ROOT="$LOCAL_BASE/tools/uv-0.11.25"
    if [[ ! -x "$UV_TOOL_ROOT/bin/uv" ]]; then
        echo "### bootstrapping uv into localscratch on the compute node"
        python3 -m venv "$UV_TOOL_ROOT"
        "$UV_TOOL_ROOT/bin/python" -m pip install --disable-pip-version-check "uv==0.11.25"
    fi
    export PATH="$UV_TOOL_ROOT/bin:$PATH"
fi
uv --version

# Installation and dependency downloads happen on this compute node. `--frozen`
# guarantees that the committed lock file, not live dependency resolution, is used.
uv sync --frozen

"$UV_PROJECT_ENVIRONMENT/bin/python" - <<'PY'
from importlib.metadata import version

for package in ("rapids-singlecell-cu12", "scanpy", "squidpy", "decoupler", "pertpy"):
    print(f"{package}=={version(package)}")
PY

set +e
"$UV_PROJECT_ENVIRONMENT/bin/python" benchmarks/comparison/run_structured.py
comparison_status=$?
set -e

# Copy only durable, compact evidence back to shared storage.
rsync -a benchmarks/comparison/results/ "$DURABLE_OUTPUT/results/"
if [[ -d benchmarks/comparison/report ]]; then
    rsync -a benchmarks/comparison/report/ "$DURABLE_OUTPUT/report/"
fi
if [[ -f benchmarks/comparison/equivalence.json ]]; then
    cp benchmarks/comparison/equivalence.json "$DURABLE_OUTPUT/equivalence.json"
fi
cp "$EQUIVALENCE_EXECUTION" "$DURABLE_OUTPUT/execution.json"
"$UV_PROJECT_ENVIRONMENT/bin/python" -m pip freeze > "$DURABLE_OUTPUT/environment.txt" 2>/dev/null || \
    uv pip freeze --python "$UV_PROJECT_ENVIRONMENT/bin/python" > "$DURABLE_OUTPUT/environment.txt"
nvidia-smi -q > "$DURABLE_OUTPUT/nvidia-smi.txt"

echo "### comparison_status=$comparison_status finished=$(date --iso-8601=seconds)"
echo "### evidence=$DURABLE_OUTPUT"
exit "$comparison_status"
