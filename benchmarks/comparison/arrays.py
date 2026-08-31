"""Store the raw CPU and GPU outputs. Comparing them is `evaluate.py`'s job.

One store per method group, one group per comparison point:

    arrays/scanpy_core_preprocessing.zarr/
        normalize_total/{reference,candidate}
        scale/{reference,candidate}

Sparse input is stored as its CSR components rather than densified — the count matrices in
this suite are ~3% dense, so densifying them would cost 40x the space for no information.
Zarr v3 with Blosc/Zstd, chunked along the cell axis, which is the axis everything here
reads along.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import zarr
from numcodecs.zarr3 import Blosc
from scipy import sparse

# Everything a run writes, unless an EQUIVALENCE_* variable overrides it.
OUT = Path(__file__).parents[2] / "out"

CELL_CHUNK = 2048


def store_path(method: str) -> Path:
    """Where this method group's arrays live. Sits beside the JSON records, not inside them."""
    root = Path(os.environ.get("EQUIVALENCE_ARRAY_DIR", OUT / "arrays"))
    return root / f"{method}.zarr"


def _write(group: zarr.Group, name: str, value) -> None:
    target = group.require_group(name)
    if sparse.issparse(value):
        csr = value.tocsr()
        target.attrs["format"] = "csr"
        target.attrs["shape"] = list(csr.shape)
        for part, data in (("data", csr.data), ("indices", csr.indices), ("indptr", csr.indptr)):
            array = np.ascontiguousarray(data)
            target.create_array(
                part,
                shape=array.shape,
                dtype=array.dtype,
                chunks=(min(len(array), CELL_CHUNK * 64),),
                compressors=[Blosc(cname="zstd", clevel=3)],
            )[:] = array
        return
    array = np.asarray(value)
    if array.dtype.kind in "OU S":
        array = array.astype(str)
    target.attrs["format"] = "dense"
    target.attrs["shape"] = list(array.shape)
    chunks = (min(array.shape[0], CELL_CHUNK), *array.shape[1:]) if array.ndim else None
    target.create_array(
        "values",
        shape=array.shape,
        dtype=array.dtype,
        chunks=chunks,
        compressors=[Blosc(cname="zstd", clevel=3)],
    )[:] = array


def capture(method: str, point: str, **arrays) -> None:
    """Record the outputs at one comparison point. Nothing is compared here.

    Arrays are keyword-only and named. By convention `reference` is the CPU side and
    `candidate` the GPU side; a comparison that needs more — trustworthiness needs the
    basis its embedding was built from — takes extra names such as `reference_basis`.
    Keyword-only because these are easy to transpose and several comparisons are
    asymmetric.
    """
    root = zarr.open_group(str(store_path(method)), mode="a")
    for name, value in arrays.items():
        if value is None:
            continue
        _write(root, f"{point}/{name}", value)


def _read(group: zarr.Group):
    if group.attrs["format"] == "csr":
        shape = tuple(group.attrs["shape"])
        return sparse.csr_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape)
    return group["values"][:]


def load_arrays(method: str, point: str, names: tuple[str, ...]):
    """Return the named arrays at one comparison point, or None if any is missing."""
    path = store_path(method)
    if not path.exists():
        return None
    root = zarr.open_group(str(path), mode="r")
    if point not in root:
        return None
    stored = root[point]
    if any(name not in stored for name in names):
        return None
    return tuple(_read(stored[name]) for name in names)


def captured_points(method: str) -> list[str]:
    """Names of every comparison point stored for this method group."""
    path = store_path(method)
    if not path.exists():
        return []
    return sorted(zarr.open_group(str(path), mode="r").group_keys())
