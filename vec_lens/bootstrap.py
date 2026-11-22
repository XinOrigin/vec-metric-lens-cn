"""Cell-level bootstrap 95% confidence intervals and floor significance checks."""

from __future__ import annotations

from typing import Any

import numpy as np


def _ci95(values: list[float]) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    return {
        "mean": round(float(np.mean(arr)), 5),
        "ci_low": round(float(np.percentile(arr, 2.5)), 5),
        "ci_high": round(float(np.percentile(arr, 97.5)), 5),
    }


def bootstrap_diagnostics(
    X_pred: np.ndarray,
    X_tgt: np.ndarray,
    X_ref: np.ndarray,
    *,
    C_pred: np.ndarray | None = None,
    C_tgt: np.ndarray | None = None,
    n_boot: int = 50,
    seed: int = 0,
) -> dict[str, Any]:
    """Compute bootstrap 95% CIs for DE correlation, centroid L2, and spatial spread error."""
    if n_boot < 5:
        raise ValueError("n_boot must be >= 5 for meaningful percentile intervals")

    rng = np.random.default_rng(seed)
    n_p, n_t, n_r = len(X_pred), len(X_tgt), len(X_ref)
    mu_ref_full = X_ref.mean(axis=0)

    expr_l2_pred: list[float] = []
    expr_l2_floor: list[float] = []
    de_r_pred: list[float] = []
    spatial_spread_err: list[float] = []

    for _ in range(n_boot):
        ip = rng.integers(0, n_p, size=n_p)
        it = rng.integers(0, n_t, size=n_t)
        ir = rng.integers(0, n_r, size=n_r)

        mu_p = X_pred[ip].mean(axis=0)
        mu_t = X_tgt[it].mean(axis=0)
        mu_r = X_ref[ir].mean(axis=0)

        expr_l2_pred.append(float(np.linalg.norm(mu_p - mu_t)))
        expr_l2_floor.append(float(np.linalg.norm(mu_r - mu_t)))

        dp = mu_p - mu_ref_full
        dt = mu_t - mu_ref_full
        if np.std(dp) > 1e-8 and np.std(dt) > 1e-8:
            de_r_pred.append(float(np.corrcoef(dp, dt)[0, 1]))
        else:
            de_r_pred.append(0.0)

        if C_pred is not None and C_tgt is not None:
            sp_p = np.std(C_pred[ip], axis=0)
            sp_t = np.std(C_tgt[it], axis=0)
            spatial_spread_err.append(float(np.linalg.norm(sp_p - sp_t)))

    pred_ci = _ci95(expr_l2_pred)
    floor_ci = _ci95(expr_l2_floor)
    # Improvement over no-change floor is significant when pred upper CI < floor lower CI
    significant_vs_floor = pred_ci["ci_high"] < floor_ci["ci_low"]

    out: dict[str, Any] = {
        "n_boot": n_boot,
        "small_sample_warning": min(n_p, n_t) < 500,
        "expr_centroid_l2": pred_ci,
        "floor_expr_centroid_l2": floor_ci,
        "beats_no_change_floor_significantly": bool(significant_vs_floor),
        "de_pearson_r": _ci95(de_r_pred),
    }
    if spatial_spread_err:
        out["spatial_spread_l2"] = _ci95(spatial_spread_err)
    return out
