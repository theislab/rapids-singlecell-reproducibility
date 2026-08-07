# CUDA comes from the pinned wheels in uv.lock, so the image needs nothing from the host
# but the driver: `docker run --gpus all`, or `apptainer run --nv`.
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/repro/.venv \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg

WORKDIR /repro
COPY pyproject.toml uv.lock ./
# Everything else in the lock file resolves to a wheel; bbknn pulls annoy, which has no
# cp312 wheel and builds from source. APT::Sandbox::User=root keeps apt from dropping to
# `_apt`, which it cannot do under a rootless builder with only one UID to map.
RUN apt-get -o APT::Sandbox::User=root update \
    && apt-get -o APT::Sandbox::User=root install -y --no-install-recommends g++ \
    && rm -rf /var/lib/apt/lists/*
# annoy's setup.py otherwise appends -march=native, which would tie the image to the CPU
# that built it. ANNOY_COMPILER_ARGS replaces the whole flag list, so drop just that one.
ENV ANNOY_COMPILER_ARGS=-D_CRT_SECURE_NO_WARNINGS,-fpermissive,-O3,-ffast-math,-fno-associative-math,-DANNOYLIB_MULTITHREADED_BUILD,-std=c++14
RUN uv sync --frozen
COPY . .

# Everything a run writes goes to /out: result records, the report, downloaded datasets and
# compilation caches. Mount it and the image itself stays read-only.
WORKDIR /out
ENV EQUIVALENCE_OUTPUT_DIR=/out/results \
    EQUIVALENCE_SUMMARY=/out/equivalence.json \
    EQUIVALENCE_EXECUTION=/out/execution.json \
    EQUIVALENCE_REPORT_DIR=/out/report \
    EQUIVALENCE_ARTIFACT_DIR=/out/report/artifacts \
    XDG_CACHE_HOME=/out/.cache \
    MPLCONFIGDIR=/out/.cache/matplotlib \
    NUMBA_CACHE_DIR=/out/.cache/numba \
    CUPY_CACHE_DIR=/out/.cache/cupy

ENTRYPOINT ["/repro/.venv/bin/python"]
CMD ["/repro/benchmarks/comparison/run_structured.py"]
