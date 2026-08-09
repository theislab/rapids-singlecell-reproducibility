"""Store the raw CPU and GPU outputs so comparisons can be made later, not just now.

A comparison script used to compute a scalar and throw the arrays away, which meant every
new question about a run — which term of `numpy.allclose` bound the failure, how the error
is distributed, whether a different metric says something else — cost another GPU run. So
the arrays are written to a Zarr store instead and the comparing happens in `evaluate.py`.

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

CELL_CHUNK = 2048


def store_path(method: str) -> Path:
    """Where this method group's arrays live. Sits beside the JSON records, not inside them."""
    root = Path(os.environ.get("EQUIVALENCE_ARRAY_DIR", Path(__file__).parent / "arrays"))
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


def capture(method: str, name: str, reference, candidate) -> None:
    """Record one comparison point: what the CPU produced and what the GPU produced.

    `reference` is the CPU side. Nothing is compared here — that is `evaluate.py`'s job.
    """
    root = zarr.open_group(str(store_path(method)), mode="a")
    _write(root, f"{name}/reference", reference)
    _write(root, f"{name}/candidate", candidate)


def _read(group: zarr.Group):
    if group.attrs["format"] == "csr":
        shape = tuple(group.attrs["shape"])
        return sparse.csr_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape)
    return group["values"][:]


def load_pair(method: str, name: str):
    """Return `(reference, candidate)` for one comparison point, or None if not captured."""
    path = store_path(method)
    if not path.exists():
        return None
    root = zarr.open_group(str(path), mode="r")
    if name not in root:
        return None
    pair = root[name]
    if "reference" not in pair or "candidate" not in pair:
        return None
    return _read(pair["reference"]), _read(pair["candidate"])


def captured_points(method: str) -> list[str]:
    """Names of every comparison point stored for this method group."""
    path = store_path(method)
    if not path.exists():
        return []
    return sorted(zarr.open_group(str(path), mode="r").group_keys())
