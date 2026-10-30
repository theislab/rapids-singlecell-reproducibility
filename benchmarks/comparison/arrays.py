from __future__ import annotations

import os
import shutil
from pathlib import Path

import numpy as np
import zarr
from numcodecs.zarr3 import Blosc
from scipy import sparse

# Everything a run writes: results/, arrays/, report/, and the run's json records.
OUT = Path(os.environ.get("EQUIVALENCE_OUT", Path(__file__).parents[2] / "out"))

CELL_CHUNK = 2048


def store_path(method: str) -> Path:
    return OUT / "arrays" / f"{method}.zarr"


def _store(group: zarr.Group, name: str, array: np.ndarray, chunks) -> None:
    group.create_array(
        name,
        shape=array.shape,
        dtype=array.dtype,
        chunks=chunks,
        compressors=[Blosc(cname="zstd", clevel=3)],
        overwrite=True,
    )[:] = array


def _write(group: zarr.Group, name: str, value) -> None:
    if type(value).__module__.startswith("cupy"):  # cupy arrays and cupyx sparse matrices
        value = value.get()
    target = group.require_group(name)
    if sparse.issparse(value):
        csr = value.tocsr()
        target.attrs["format"] = "csr"
        target.attrs["shape"] = list(csr.shape)
        for part in ("data", "indices", "indptr"):
            array = np.ascontiguousarray(getattr(csr, part))
            _store(target, part, array, (max(1, min(len(array), CELL_CHUNK * 64)),))
        return
    array = np.asarray(value)
    if array.dtype.kind in "OUS":
        array = array.astype(str)
    target.attrs["format"] = "dense"
    target.attrs["shape"] = list(array.shape)
    _store(target, "values", array, (max(1, min(array.shape[0], CELL_CHUNK)), *array.shape[1:]))


_STARTED: set[str] = set()


def capture(method: str, point: str, **arrays) -> None:
    """Store the outputs at one comparison point; `reference` is CPU, `candidate` GPU."""
    if method not in _STARTED:  # each run of a script starts from an empty store
        shutil.rmtree(store_path(method), ignore_errors=True)
        _STARTED.add(method)
    root = zarr.open_group(str(store_path(method)), mode="a")
    point = point.replace("/", "_")  # "/" would nest the point in zarr
    for name, value in arrays.items():
        _write(root, f"{point}/{name}", value)


def _read(group: zarr.Group):
    if group.attrs["format"] == "csr":
        shape = tuple(group.attrs["shape"])
        return sparse.csr_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape)
    return group["values"][:]


def load_arrays(method: str, point: str, names: tuple[str, ...]):
    """Return the named arrays at one comparison point, or None if any is missing."""
    stored = zarr.open_group(str(store_path(method)), mode="r")[point]
    if any(name not in stored for name in names):
        return None
    return tuple(_read(stored[name]) for name in names)


def captured_points(method: str) -> list[str]:
    """Names of every comparison point stored for this method group."""
    path = store_path(method)
    if not path.exists():
        return []
    return sorted(zarr.open_group(str(path), mode="r").group_keys())
