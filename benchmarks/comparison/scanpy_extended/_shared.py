from __future__ import annotations

import scanpy as sc


def pbmc68k():
    return sc.datasets.pbmc68k_reduced()
