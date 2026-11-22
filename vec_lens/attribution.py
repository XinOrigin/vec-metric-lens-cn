"""Per-gene and leave-one-celltype-out metric attribution for VEC predictions."""

from __future__ import annotations

from typing import Any

import numpy as np


def explain_gene_de(
    X_pred: np.ndarray,
    X_tgt: np.ndarray,
    X_ref: np.ndarray,
    genes: list[str],
    *,
    top_k: int = 10,
) -> dict[str, Any]:
    """Identify genes driving differential-expression (DE) magnitude & direction errors."""
    mu_pred = X_pred.mean(axis=0)
    mu_tgt = X_tgt.mean(axis=0)
    mu_ref = X_ref.mean(axis=0)

    pred_delta = mu_pred - mu_ref
    tgt_delta = mu_tgt - mu_ref
    abs_err = np.abs(pred_delta - tgt_delta)

    # Pearson correlation of DE deltas (aligned with de_score/de_direction intuition)
    if np.std(pred_delta) > 1e-8 and np.std(tgt_delta) > 1e-8:
        de_corr = float(np.corrcoef(pred_delta, tgt_delta)[0, 1])
    else:
        de_corr = 0.0

    nonzero_mask = np.abs(tgt_delta) > 1e-4
    sign_flipped = nonzero_mask & ((pred_delta * tgt_delta) < 0)
    sign_agree_rate = (
        float(np.mean((pred_delta[nonzero_mask] * tgt_delta[nonzero_mask]) > 0))
        if np.any(nonzero_mask)
        else 1.0
    )

    err_order = np.argsort(-abs_err)[:top_k]
    top_error_genes = [
        {
            "gene": genes[int(i)],
            "pred_delta": round(float(pred_delta[i]), 5),
            "target_delta": round(float(tgt_delta[i]), 5),
            "abs_error": round(float(abs_err[i]), 5),
            "direction_flipped": bool(sign_flipped[i]),
        }
        for i in err_order
    ]

    flip_indices = np.where(sign_flipped)[0]
    flip_sorted = flip_indices[np.argsort(-abs_err[flip_indices])][:top_k]
    top_flipped_genes = [
        {
            "gene": genes[int(i)],
            "pred_delta": round(float(pred_delta[i]), 5),
            "target_delta": round(float(tgt_delta[i]), 5),
            "abs_error": round(float(abs_err[i]), 5),
        }
        for i in flip_sorted
    ]

    return {
        "de_pearson_r": round(de_corr, 4),
        "de_sign_agreement": round(sign_agree_rate, 4),
        "n_flipped_genes": int(np.sum(sign_flipped)),
        "mean_abs_de_error": round(float(np.mean(abs_err)), 5),
        "top_error_genes": top_error_genes,
        "top_flipped_genes": top_flipped_genes,
    }


def _centroid_mismatch(
    X_a: np.ndarray,
    X_b: np.ndarray,
    C_a: np.ndarray | None,
    C_b: np.ndarray | None,
) -> tuple[float, float | None]:
    expr_l2 = float(np.linalg.norm(X_a.mean(axis=0) - X_b.mean(axis=0)))
    spatial_l2: float | None = None
    if C_a is not None and C_b is not None and len(C_a) > 0 and len(C_b) > 0:
        spatial_l2 = float(np.linalg.norm(C_a.mean(axis=0) - C_b.mean(axis=0)))
    return expr_l2, spatial_l2


def explain_celltype_logo(
    X_pred: np.ndarray,
    X_tgt: np.ndarray,
    ct_pred: np.ndarray,
    ct_tgt: np.ndarray,
    *,
    C_pred: np.ndarray | None = None,
    C_tgt: np.ndarray | None = None,
    top_k: int = 8,
) -> dict[str, Any]:
    """Leave-one-celltype-out attribution on expression & 3D spatial centroid shift."""
    base_expr, base_spatial = _centroid_mismatch(X_pred, X_tgt, C_pred, C_tgt)
    groups = sorted(set(ct_pred.tolist()) | set(ct_tgt.tolist()))
    impacts: list[dict[str, Any]] = []

    for grp in groups:
        keep_p = ct_pred != grp
        keep_t = ct_tgt != grp
        if not np.any(keep_p) or not np.any(keep_t):
            continue
        Cp_sub = C_pred[keep_p] if C_pred is not None else None
        Ct_sub = C_tgt[keep_t] if C_tgt is not None else None
        expr_wo, spatial_wo = _centroid_mismatch(X_pred[keep_p], X_tgt[keep_t], Cp_sub, Ct_sub)
        # Positive delta means removing this cell type REDUCES mismatch -> major error contributor
        expr_contrib = base_expr - expr_wo
        spatial_contrib = (
            round(base_spatial - spatial_wo, 5)
            if (base_spatial is not None and spatial_wo is not None)
            else None
        )
        impacts.append(
            {
                "celltype": grp,
                "n_pred": int(np.sum(~keep_p)),
                "n_target": int(np.sum(~keep_t)),
                "expr_l2_reduction_when_removed": round(expr_contrib, 5),
                "spatial_l2_reduction_when_removed": spatial_contrib,
            }
        )

    impacts.sort(key=lambda r: -abs(r["expr_l2_reduction_when_removed"]))
    return {
        "baseline_expr_centroid_l2": round(base_expr, 5),
        "baseline_spatial_centroid_l2": round(base_spatial, 5) if base_spatial is not None else None,
        "celltypes": impacts[:top_k],
    }
