# vec-metric-lens-cn

面向 **Virtual Embryo Challenge (`veckit`)** 的中英双语指标归因与生物学合理性诊断工具。

## 为什么需要 `vec-metric-lens-cn`？

- `vec-preflight-cn` 负责**提交前格式与泄漏门禁**；
- `vec-baseline-bench-cn` 负责**无泄漏基准评测与 3D 对比**；
- **`vec-metric-lens-cn` 负责解释分数为什么变化**：定位误差最大的基因、识别方向反转基因（`direction_flipped`）、执行细胞类型留一归因（Leave-One-Group-Out）、输出 95% Bootstrap 置信区间与小样本告警，并在提交前验证 Marker 基因生物学方向一致性。

## 快速安装 / Quick Start

```bash
pip install -e .
# 若需同步输出 veckit 官方面板分数：
pip install git+https://github.com/aristoteleo/vec_baselines.git
```

国内镜像加速：

```bash
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple -e .
```

## 使用示例 / Usage

```bash
vec-lens explain \
  --task T3 \
  --input pred.h5ad \
  --target mab21l2_ko.h5ad \
  --reference wt.h5ad \
  --obs-key celltype \
  --n-boot 50 \
  --with-veckit \
  --veckit-dir ../veckit \
  --out-dir lens_out/
```

详见中文操作指南：[docs/LENS_ZH.md](docs/LENS_ZH.md)。
