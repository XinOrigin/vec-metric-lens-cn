# vec-metric-lens-cn 中文诊断指南

## 0. 环境安装

推荐使用 Python `>= 3.10`：

```bash
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple -e .
pip install git+https://github.com/aristoteleo/vec_baselines.git
```

## 1. 输出内容说明

运行 `vec-lens explain --out-dir lens_out/` 后会生成：

- `lens_report.json`：包含 `gene_de`（Top 误差基因与反转基因）、`celltype_logo`（细胞类型留一法对表达/空间质心 L2 误差贡献）、`bootstrap`（95% 置信区间及是否显著优于无变化基线 `beats_no_change_floor_significantly`）、`biology`（Marker 基因上下调一致性、细胞比例 TVD、稀疏度差异）以及可选的 `veckit` 原始指标。
- `lens_report.html`：完全离线、零 CDN 依赖的可视化中文报告。

> **免责声明**：本工具仅读取用户本地传入的 `.h5ad` 文件，不访问任何隐藏测试集；离线诊断分数仅作本地调优参考。
