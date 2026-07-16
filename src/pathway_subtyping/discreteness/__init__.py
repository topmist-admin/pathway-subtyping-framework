"""
Discreteness-aware subtype validation (added in v0.8.0).

This subpackage corrects and extends the stability gate. The original
``ValidationGates`` bootstrap-stability null permuted each pathway column
independently, which tests the hypothesis *"the pathways are mutually
independent"* — **not** *"there are discrete clusters."* Those diverge exactly at
small n: a single correlated Gaussian blob or a 1-D continuous gradient (tumor
purity, immune infiltration, proliferation) has no discrete clusters, yet a
mixture model reproducibly bisects it every bootstrap (high observed ARI) while
the feature-permuted null collapses (low ARI), so a continuum is falsely
certified a "reproducible subtype."

``DiscretenessGateA`` keeps the identical observed statistic (mean bootstrap-ARI
at the fixed selected k) and replaces the reference with a **discreteness** null:

  * primary  — single-Gaussian (SigClust) reference: is the partition more stable
    than partitions of data drawn from one Gaussian of the same covariance?
    (Liu, Hayes, Nobel & Marron, JASA 2008, doi:10.1198/016214508000000454)
  * complementary — gap statistic (Tibshirani, Walther & Hastie, JRSS B 2001,
    doi:10.1111/1467-9868.00293)
  * corroborating — Hartigan's dip test of unimodality (Ann. Statist. 1985,
    doi:10.1214/aos/1176346577); optional, requires the ``diptest`` extra.

The old independence null is **retained but demoted** to a confound / marginal
control (reported, no longer the gate). Small-n hardening (PCA to d ≪ n, fixed k
across observed and reference, silhouette-based k-stability routing) is applied
identically to observed data and every reference draw.

``ReframedMembershipGate`` reframes the conformal membership gate as **assignment
sharpness conditional on a partition already certified by Gates A and B** — not a
standalone guarantee — at a single 0.90 operating point with a coverage CI and a
minimum-n rule.

Import path::

    from pathway_subtyping.discreteness import DiscretenessGateA, ReframedMembershipGate
"""
from .gate_a_discreteness_null import (
    DiscretenessGateA,
    GateAv2Result,
    gap_statistic,
    reduced_dim,
    reduce_scores,
)
from .gate_c_conformal_membership import ConformalMembershipGate
from .gate_c_reframed import ReframedMembershipGate, GateCReframedResult

__all__ = [
    "DiscretenessGateA",
    "GateAv2Result",
    "ReframedMembershipGate",
    "GateCReframedResult",
    "ConformalMembershipGate",
    "gap_statistic",
    "reduced_dim",
    "reduce_scores",
]
