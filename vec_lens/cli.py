"""CLI entrypoint for vec-metric-lens-cn (`vec-lens explain`)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from .attribution import explain_celltype_logo, explain_gene_de
from .biology import check_biology_sanity
from .bootstrap import bootstrap_diagnostics
from .report import render_html_report
from .scoring import load_h5ad_bundle, run_veckit_score


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="vec-lens", description="Bilingual VEC score attribution & diagnostics")
    sub = p.add_subparsers(dest="command", required=True)

    ex = sub.add_parser("explain", help="Explain metric drivers, bootstrap CIs, and biological sanity")
    ex.add_argument("--task", required=True, choices=["T1", "T2", "T3", "t1", "t2", "t3"])
    ex.add_argument("--input", required=True, type=Path, help="Predicted .h5ad path")
    ex.add_argument("--target", required=True, type=Path, help="Local pseudo-target .h5ad path")
    ex.add_argument("--reference", type=Path, default=None, help="Reference/WT .h5ad path (defaults to --input)")
    ex.add_argument("--obs-key", default="celltype", help="obs column for cell-type attribution")
    ex.add_argument("--markers", type=Path, default=None, help="JSON file with {'up': [...], 'down': [...]}")
    ex.add_argument("--n-boot", type=int, default=50, help="Bootstrap resamples (default: 50)")
    ex.add_argument("--seed", type=int, default=0)
    ex.add_argument("--allow-reorder", action="store_true")
    ex.add_argument("--setting", default="heart")
    ex.add_argument("--with-veckit", action="store_true", help="Also run veckit.score() panel")
    ex.add_argument("--veckit-dir", type=Path, default=None, help="Path to local veckit checkout")
    ex.add_argument("--out-dir", type=Path, default=None, help="Directory to write lens_report.json & .html")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    if args.command != "explain":
        return 2

    task = args.task.upper()
    bundle = load_h5ad_bundle(
        args.input,
        args.target,
        args.reference,
        obs_key=args.obs_key,
        allow_reorder=args.allow_reorder,
    )

    payload = {
        "task": task,
        "gene_de": explain_gene_de(bundle.X_pred, bundle.X_tgt, bundle.X_ref, bundle.genes),
        "celltype_logo": explain_celltype_logo(
            bundle.X_pred,
            bundle.X_tgt,
            bundle.ct_pred,
            bundle.ct_tgt,
            C_pred=bundle.C_pred,
            C_tgt=bundle.C_tgt,
        ),
        "bootstrap": bootstrap_diagnostics(
            bundle.X_pred,
            bundle.X_tgt,
            bundle.X_ref,
            C_pred=bundle.C_pred,
            C_tgt=bundle.C_tgt,
            n_boot=args.n_boot,
            seed=args.seed,
        ),
        "biology": check_biology_sanity(
            bundle.X_pred,
            bundle.X_tgt,
            bundle.X_ref,
            bundle.genes,
            bundle.ct_pred,
            bundle.ct_tgt,
            celltype_defaulted=bundle.celltype_defaulted,
            markers_path=args.markers,
        ),
        "veckit": (
            run_veckit_score(
                task,
                args.input,
                args.target,
                reference_path=args.reference,
                setting=args.setting,
                seed=args.seed,
                allow_reorder=args.allow_reorder,
                veckit_dir=args.veckit_dir,
            )
            if args.with_veckit
            else None
        ),
    }

    if args.out_dir is not None:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        (args.out_dir / "lens_report.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        (args.out_dir / "lens_report.html").write_text(
            render_html_report(payload), encoding="utf-8"
        )
    else:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
