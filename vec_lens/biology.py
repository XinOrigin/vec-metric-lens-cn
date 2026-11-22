"""Biological sanity checks: marker gene directionality, proportions, and sparsity."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


def _tvd(p: dict[str, float], q: dict[str, float]) -> float:
    keys = set(p) | set(q)
    return 0.5 * sum(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in keys)


def _proportions(labels: np.ndarray) -> dict[str, float]:
    n = max(len(labels), 1)
    vals, counts = np.unique(labels, return_counts=True)
    return {str(v): float(c) / n for v, c in zip(vals, counts)}


def check_biology_sanity(
    X_pred: np.ndarray,
    X_tgt: np.ndarray,
    X_ref: np.ndarray,
    genes: list[str],
    ct_pred: np.ndarray,
    ct_tgt: np.ndarray,
    *,
    celltype_defaulted: bool = False,
    markers_path: Path | None = None,
) -> dict[str, Any]:
    """Verify marker directionality, cell-type proportion shift, and sparsity parity."""
    gene_idx = {g: idx for idx, g in enumerate(genes)}

    pred_zero_rate = float(np.mean(X_pred == 0))
    tgt_zero_rate = float(np.mean(X_tgt == 0))

    p_prop = _proportions(ct_pred)
    t_prop = _proportions(ct_tgt)
    prop_tvd = None if celltype_defaulted else round(_tvd(p_prop, t_prop), 4)

    marker_checks: list[dict[str, Any]] = []
    if markers_path is not None:
        spec = json.loads(Path(markers_path).read_text(encoding="utf-8"))
        mu_pred = X_pred.mean(axis=0)
        mu_tgt = X_tgt.mean(axis=0)
        mu_ref = X_ref.mean(axis=0)
        for expected_dir, names in [("up", spec.get("up", [])), ("down", spec.get("down", []))]:
            for g in names:
                if g not in gene_idx:
                    marker_checks.append({"gene": g, "expected": expected_dir, "status": "missing_in_panel"})
                    continue
                i = gene_idx[g]
                p_delta = float(mu_pred[i] - mu_ref[i])
                t_delta = float(mu_tgt[i] - mu_ref[i])
                ok = (p_delta > 0) if expected_dir == "up" else (p_delta < 0)
                marker_checks.append(
                    {
                        "gene": g,
                        "expected": expected_dir,
                        "pred_delta": round(p_delta, 5),
                        "target_delta": round(t_delta, 5),
                        "ok": bool(ok),
                    }
                )

    return {
        "celltype_defaulted": celltype_defaulted,
        "celltype_proportion_tvd": prop_tvd,
        "pred_sparsity": round(pred_zero_rate, 4),
        "target_sparsity": round(tgt_zero_rate, 4),
        "sparsity_abs_diff": round(abs(pred_zero_rate - tgt_zero_rate), 4),
        "marker_checks": marker_checks,
    }
