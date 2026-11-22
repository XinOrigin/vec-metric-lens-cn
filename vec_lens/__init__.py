"""vec-metric-lens-cn: bilingual diagnostic & attribution lens for Virtual Embryo Challenge."""

from .attribution import explain_celltype_logo, explain_gene_de
from .biology import check_biology_sanity
from .bootstrap import bootstrap_diagnostics
from .scoring import load_h5ad_bundle, run_veckit_score

__all__ = [
    "bootstrap_diagnostics",
    "check_biology_sanity",
    "explain_celltype_logo",
    "explain_gene_de",
    "load_h5ad_bundle",
    "run_veckit_score",
]
