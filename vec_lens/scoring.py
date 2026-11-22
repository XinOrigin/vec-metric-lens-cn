"""Clean bridge to veckit offline scoring and validated AnnData loading."""

from __future__ import annotations

import importlib
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse


@dataclass
class BundleData:
    X_pred: np.ndarray
    X_tgt: np.ndarray
    X_ref: np.ndarray
    C_pred: np.ndarray | None
    C_tgt: np.ndarray | None
    genes: list[str]
    ct_pred: np.ndarray
    ct_tgt: np.ndarray
    celltype_defaulted: bool


def _to_dense_float32(X: Any) -> np.ndarray:
    arr = X.toarray() if sparse.issparse(X) else np.asarray(X)
    return np.asarray(arr, dtype=np.float32)


def _extract_coords(adata: ad.AnnData) -> np.ndarray | None:
    if "spatial_3D" not in adata.obsm:
        return None
    C = np.asarray(adata.obsm["spatial_3D"], dtype=np.float32)
    if C.ndim != 2 or C.shape[1] < 3:
        return None
    return C[:, :3]


def load_h5ad_bundle(
    input_path: Path,
    target_path: Path,
    reference_path: Path | None = None,
    *,
    obs_key: str = "celltype",
    allow_reorder: bool = False,
) -> BundleData:
    """Load prediction, target, and reference/WT AnnData files cleanly."""
    a_pred = ad.read_h5ad(input_path)
    a_tgt = ad.read_h5ad(target_path)
    a_ref = ad.read_h5ad(reference_path) if reference_path is not None else a_pred

    genes = [str(g) for g in a_pred.var_names]
    for label, obj in [("target", a_tgt), ("reference", a_ref)]:
        obj_genes = [str(g) for g in obj.var_names]
        if obj_genes != genes:
            if allow_reorder and set(obj_genes) == set(genes) and len(obj_genes) == len(genes):
                if label == "target":
                    a_tgt = a_tgt[:, genes].copy()
                else:
                    a_ref = a_ref[:, genes].copy()
            else:
                raise ValueError(f"{label} var_names do not match prediction gene panel")

    ct_defaulted = False
    if obs_key in a_pred.obs and obs_key in a_tgt.obs:
        ct_pred = np.asarray(a_pred.obs[obs_key]).astype(str)
        ct_tgt = np.asarray(a_tgt.obs[obs_key]).astype(str)
    else:
        ct_defaulted = True
        ct_pred = np.array(["NA"] * a_pred.n_obs)
        ct_tgt = np.array(["NA"] * a_tgt.n_obs)

    return BundleData(
        X_pred=_to_dense_float32(a_pred.X),
        X_tgt=_to_dense_float32(a_tgt.X),
        X_ref=_to_dense_float32(a_ref.X),
        C_pred=_extract_coords(a_pred),
        C_tgt=_extract_coords(a_tgt),
        genes=genes,
        ct_pred=ct_pred,
        ct_tgt=ct_tgt,
        celltype_defaulted=ct_defaulted,
    )


def run_veckit_score(
    task: str,
    input_path: Path,
    target_path: Path,
    *,
    reference_path: Path | None = None,
    wt_path: Path | None = None,
    setting: str = "heart",
    seed: int = 0,
    allow_reorder: bool = False,
    veckit_dir: Path | None = None,
) -> dict[str, Any]:
    """Call veckit.score() without duplicating official metric formulas."""
    if veckit_dir is not None:
        vdir = str(Path(veckit_dir).resolve())
        if vdir not in sys.path:
            sys.path.insert(0, vdir)
    try:
        veckit_mod = importlib.import_module("veckit")
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "veckit is not installed. Install via "
            "`pip install git+https://github.com/aristoteleo/vec_baselines.git` "
            "or pass --veckit-dir PATH."
        ) from exc

    kw: dict[str, Any] = {
        "seed": seed,
        "allow_reorder": allow_reorder,
        "setting": setting,
    }
    if task.upper() == "T3":
        kw["wt"] = wt_path or reference_path
    else:
        kw["reference"] = reference_path

    return veckit_mod.score(task=task.upper(), input=input_path, target=target_path, **kw)
