"""Offline self-contained bilingual (CN/EN) HTML report generator (no CDN)."""

from __future__ import annotations

import html
import json
from typing import Any


def render_html_report(payload: dict[str, Any]) -> str:
    """Render a standalone offline HTML diagnostic report."""
    task = html.escape(str(payload.get("task", "")))
    de = payload.get("gene_de", {})
    boot = payload.get("bootstrap", {})
    bio = payload.get("biology", {})
    logo = payload.get("celltype_logo", {})
    veckit_res = payload.get("veckit")

    warnings: list[str] = []
    if bio.get("celltype_defaulted"):
        warnings.append("观测到 obs['celltype'] 缺失，已自动回退为 'NA'，细胞类型归因仅供参考。")
    if boot.get("small_sample_warning"):
        warnings.append("样本细胞数 < 500：小样本下指标波动较大，请注意 95% 置信区间宽度。")
    if not boot.get("beats_no_change_floor_significantly", False):
        warnings.append("Δ vs floor not significant（95% CI 与无变化基线存在重叠）。")
    warnings.append("注意：离线本地诊断分数不代表隐藏测试集最终成绩（仅使用用户本地传入文件）。")

    warn_html = "".join(f"<li>{html.escape(w)}</li>" for w in warnings)

    gene_rows = "".join(
        f"<tr><td>{html.escape(g['gene'])}</td><td>{g['pred_delta']:+.4f}</td>"
        f"<td>{g['target_delta']:+.4f}</td><td>{g['abs_error']:.4f}</td>"
        f"<td>{'FLIP' if g['direction_flipped'] else 'OK'}</td></tr>"
        for g in de.get("top_error_genes", [])
    )

    ct_rows = "".join(
        f"<tr><td>{html.escape(c['celltype'])}</td><td>{c['n_pred']}</td><td>{c['n_target']}</td>"
        f"<td>{c['expr_l2_reduction_when_removed']:+.5f}</td>"
        f"<td>{c['spatial_l2_reduction_when_removed'] if c['spatial_l2_reduction_when_removed'] is not None else '—'}</td></tr>"
        for c in logo.get("celltypes", [])
    )

    veckit_pre = (
        f"<pre>{html.escape(json.dumps(veckit_res, ensure_ascii=False, indent=2))}</pre>"
        if veckit_res is not None
        else "<p>未启用 veckit 官方面板评分（可使用 <code>--with-veckit</code> 开启）。</p>"
    )

    return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>vec-metric-lens-cn 诊断报告 ({task})</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 2rem auto; max-width: 960px; color: #1f2937; line-height: 1.5; }}
h1, h2 {{ border-bottom: 1px solid #e5e7eb; padding-bottom: 0.3rem; }}
.warn {{ background: #fffbeb; border: 1px solid #f59e0b; padding: 0.8rem 1rem; border-radius: 6px; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 0.8rem; margin: 1rem 0; }}
.card {{ background: #f9fafb; border: 1px solid #e5e7eb; border-radius: 6px; padding: 0.8rem; }}
table {{ width: 100%; border-collapse: collapse; margin: 1rem 0; font-size: 0.9rem; }}
th, td {{ border: 1px solid #e5e7eb; padding: 0.45rem 0.6rem; text-align: left; }}
th {{ background: #f3f4f6; }}
pre {{ background: #111827; color: #f3f4f6; padding: 0.8rem; border-radius: 6px; overflow-x: auto; font-size: 0.82rem; }}
</style>
</head>
<body>
<h1>vec-metric-lens-cn 诊断报告 / Diagnostic Report ({task})</h1>
<div class="warn"><strong>注意事项 / Caveats</strong><ul>{warn_html}</ul></div>
<div class="grid">
  <div class="card"><strong>DE Pearson r</strong><div>{de.get('de_pearson_r')} (95% CI: {boot.get('de_pearson_r', {}).get('ci_low')} ~ {boot.get('de_pearson_r', {}).get('ci_high')})</div></div>
  <div class="card"><strong>Sign Agreement</strong><div>{de.get('de_sign_agreement')} (Flipped: {de.get('n_flipped_genes')})</div></div>
  <div class="card"><strong>Expr Centroid L2</strong><div>{boot.get('expr_centroid_l2', {}).get('mean')} vs Floor {boot.get('floor_expr_centroid_l2', {}).get('mean')}</div></div>
  <div class="card"><strong>Sparsity (Pred / Target)</strong><div>{bio.get('pred_sparsity')} / {bio.get('target_sparsity')}</div></div>
</div>
<h2>1. Top 差异表达误差基因 / Top DE Error Genes</h2>
<table><thead><tr><th>Gene</th><th>Pred Δ</th><th>Target Δ</th><th>Abs Err</th><th>Dir</th></tr></thead><tbody>{gene_rows}</tbody></table>
<h2>2. 细胞类型留一归因 / Leave-One-Group-Out Attribution</h2>
<table><thead><tr><th>Cell Type</th><th>n_pred</th><th>n_target</th><th>Expr L2 Δ</th><th>Spatial L2 Δ</th></tr></thead><tbody>{ct_rows}</tbody></table>
<h2>3. veckit 原始指标 / veckit Raw Metrics</h2>
{veckit_pre}
</body>
</html>
"""
