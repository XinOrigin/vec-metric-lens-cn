"""Unit and integration tests for vec-metric-lens-cn."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from vec_lens.cli import main as cli_main

VECKIT_DIR = Path(__file__).resolve().parents[2] / "veckit"


def _make_adata(X: np.ndarray, genes: list[str], celltypes: list[str] | None = None, coords: np.ndarray | None = None) -> ad.AnnData:
    obs = {"celltype": celltypes} if celltypes is not None else {}
    a = ad.AnnData(X=X.astype(np.float32), obs=obs)
    a.var_names = genes
    if coords is not None:
        a.obsm["spatial_3D"] = coords.astype(np.float32)
    return a


def test_explain_end_to_end_synthetic_with_markers_and_html(tmp_path: Path):
    rng = np.random.default_rng(42)
    genes = ["Mab21l2", "Gata4", "Nkx2-5", "Tbx5"]
    cts = ["Neural"] * 20 + ["Cardiac"] * 20

    X_ref = rng.uniform(1.0, 2.0, size=(40, 4))
    X_tgt = X_ref.copy()
    X_tgt[:, 0] -= 0.8  # Mab21l2 down in KO target
    X_tgt[:, 1] += 0.9  # Gata4 up

    # Prediction matches Mab21l2 down and Gata4 up closely, flips Nkx2-5 slightly
    X_pred = X_tgt + rng.normal(0.0, 0.05, size=(40, 4))
    C_ref = rng.normal(0.0, 1.0, size=(40, 3))
    C_tgt = C_ref + 0.2

    p_ref = tmp_path / "ref.h5ad"
    p_tgt = tmp_path / "tgt.h5ad"
    p_pred = tmp_path / "pred.h5ad"
    _make_adata(X_ref, genes, cts, C_ref).write_h5ad(p_ref)
    _make_adata(X_tgt, genes, cts, C_tgt).write_h5ad(p_tgt)
    _make_adata(X_pred, genes, cts, C_tgt).write_h5ad(p_pred)

    markers = tmp_path / "markers.json"
    markers.write_text(json.dumps({"down": ["Mab21l2"], "up": ["Gata4"]}), encoding="utf-8")

    out_dir = tmp_path / "lens_out"
    rc = cli_main([
        "explain",
        "--task", "T3",
        "--input", str(p_pred),
        "--target", str(p_tgt),
        "--reference", str(p_ref),
        "--markers", str(markers),
        "--n-boot", "20",
        "--out-dir", str(out_dir),
    ])
    assert rc == 0
    rep = json.loads((out_dir / "lens_report.json").read_text(encoding="utf-8"))
    assert rep["gene_de"]["de_pearson_r"] > 0.9
    assert rep["bootstrap"]["beats_no_change_floor_significantly"] is True
    assert rep["biology"]["celltype_defaulted"] is False
    assert all(m["ok"] for m in rep["biology"]["marker_checks"])
    html_text = (out_dir / "lens_report.html").read_text(encoding="utf-8")
    assert "vec-metric-lens-cn" in html_text
    assert "Mab21l2" in html_text


def test_explain_flags_missing_celltype_column(tmp_path: Path):
    genes = ["G1", "G2", "G3"]
    X = np.ones((10, 3), dtype=np.float32)
    p_pred = tmp_path / "p.h5ad"
    p_tgt = tmp_path / "t.h5ad"
    _make_adata(X, genes, celltypes=None).write_h5ad(p_pred)
    _make_adata(X + 0.5, genes, celltypes=None).write_h5ad(p_tgt)

    out_dir = tmp_path / "out_no_ct"
    rc = cli_main([
        "explain", "--task", "T1",
        "--input", str(p_pred), "--target", str(p_tgt),
        "--n-boot", "10", "--out-dir", str(out_dir),
    ])
    assert rc == 0
    rep = json.loads((out_dir / "lens_report.json").read_text(encoding="utf-8"))
    assert rep["biology"]["celltype_defaulted"] is True
    assert rep["biology"]["celltype_proportion_tvd"] is None


@pytest.mark.skipif(not (VECKIT_DIR / "data" / "sample_8.5.h5ad").exists(), reason="local veckit samples absent")
def test_explain_with_real_veckit_samples(tmp_path: Path):
    data = VECKIT_DIR / "data"
    out_dir = tmp_path / "real_out"
    rc = cli_main([
        "explain",
        "--task", "T1",
        "--input", str(data / "sample_9.5.h5ad"),
        "--target", str(data / "sample_9.5.h5ad"),
        "--reference", str(data / "sample_8.5.h5ad"),
        "--n-boot", "10",
        "--with-veckit",
        "--veckit-dir", str(VECKIT_DIR),
        "--out-dir", str(out_dir),
    ])
    assert rc == 0
    rep = json.loads((out_dir / "lens_report.json").read_text(encoding="utf-8"))
    assert rep["gene_de"]["de_pearson_r"] == pytest.approx(1.0, abs=1e-4)
    assert rep["veckit"]["metrics"]["de_score"] > 0.4
    assert "lens_report.html" in [p.name for p in out_dir.iterdir()]
