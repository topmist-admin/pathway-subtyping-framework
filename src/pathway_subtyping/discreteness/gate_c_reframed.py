"""
Gate C, reframed (v2).

WHAT CHANGED AND WHY
--------------------
The original Gate C reported split-conformal coverage of the GMM's own hard
cluster label at three targets (0.80/0.90/0.95) and a per-patient "confident
fraction". Two problems the reviewer raised (brief section 4):

  1. CIRCULARITY OF THE TARGET. Conformal prediction guarantees coverage of
     whatever label you feed it. Here the label is the GMM assignment, so
     "coverage of the assignment" is guaranteed by construction whenever the
     conformal machinery is calibrated -- it certifies the *procedure*, not that
     the subtype is real. Gate C therefore cannot, on its own, tell a real
     subtype from a sliced continuum. It must be read as ASSIGNMENT SHARPNESS
     **conditional on** a partition that Gates A (discreteness) and B
     (replication) have already certified -- not as a standalone guarantee.

  2. THREE MOVING TARGETS. Reporting 0.80/0.90/0.95 invited target-shopping.
     We fix a SINGLE operating point (0.90) agreed in advance, and report the
     empirically achieved coverage with a confidence interval across splits,
     plus the sharpness (confident/singleton fraction) at that one target.

This module wraps the existing ConformalMembershipGate, evaluates only at the
0.90 target, forms a CI on achieved coverage from the per-split values, applies
a MIN-N RULE (below a floor the 3-way split cannot be trusted -> reported as
under-powered rather than pass/fail), and labels the whole result explicitly as
CONDITIONAL on Gate A and Gate B.

Research use only.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from .gate_c_conformal_membership import ConformalMembershipGate

# Minimum n for a trustworthy 3-way split-conformal partition:
# train (>=k*2, ~45%) + calibration (>=10) + test (>=5). With cal>=10 fixed,
# n must clear ~30 for train and test folds to both be usable; we set the floor
# at 30 and additionally require the achieved calibration fold to have held.
MIN_N_GATE_C = 30
PRIMARY_COVERAGE = 0.90
COVERAGE_TOL = 0.05  # |achieved - target| <= tol  => "coverage achieved"


@dataclass
class GateCReframedResult:
    tumor: str
    n: int
    k: int
    conditional_on: str
    target_coverage: float
    achieved_coverage: float
    coverage_ci95: List[float]
    coverage_achieved: bool  # achieved within tol of the single target
    sharpness_confident_fraction: float  # singleton-set fraction at 0.90
    sharpness_ci95: List[float]
    mean_set_size: float
    n_splits_used: int
    under_powered: bool  # n < MIN_N_GATE_C
    interpretation: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def _boot_ci(vals: np.ndarray, seed: int = 42, nb: int = 2000) -> List[float]:
    vals = np.asarray(vals, float)
    if len(vals) < 2:
        return [float("nan"), float("nan")]
    rng = np.random.default_rng(seed)
    means = [rng.choice(vals, len(vals), replace=True).mean() for _ in range(nb)]
    return [round(float(np.percentile(means, 2.5)), 4), round(float(np.percentile(means, 97.5)), 4)]


class ReframedMembershipGate:
    """Gate C at a single 0.90 operating point, with a coverage CI, a min-n
    rule, and explicit conditional-on-A/B framing."""

    def __init__(
        self,
        seed: int = 42,
        n_splits: int = 80,
        primary_coverage: float = PRIMARY_COVERAGE,
        tol: float = COVERAGE_TOL,
        min_n: int = MIN_N_GATE_C,
    ):
        self.seed = seed
        self.n_splits = n_splits
        self.primary_coverage = primary_coverage
        self.tol = tol
        self.min_n = min_n

    def run(
        self,
        tumor: str,
        pathway_scores: pd.DataFrame,
        cluster_labels: np.ndarray,
        k: int,
        gate_a_pass: Optional[bool] = None,
        gate_b_pass: Optional[bool] = None,
    ) -> GateCReframedResult:
        n = pathway_scores.shape[0]
        base = ConformalMembershipGate(
            seed=self.seed,
            n_splits=self.n_splits,
            target_coverages=(self.primary_coverage,),
            primary_coverage=self.primary_coverage,
            calib_tol=self.tol,
        )
        res = base.membership_gate(tumor, pathway_scores, cluster_labels, k)
        key = f"{self.primary_coverage:.2f}"
        cov_splits = np.asarray(base._last_cov_per_split[key], float)
        conf_splits = np.asarray(base._last_conf_per_split[key], float)

        achieved = float(np.mean(cov_splits)) if len(cov_splits) else float("nan")
        cov_ci = _boot_ci(cov_splits, self.seed)
        sharp = float(np.mean(conf_splits)) if len(conf_splits) else float("nan")
        sharp_ci = _boot_ci(conf_splits, self.seed)
        cov_ok = bool(abs(achieved - self.primary_coverage) <= self.tol)
        under = bool(n < self.min_n)

        cond_bits = []
        if gate_a_pass is not None:
            cond_bits.append(f"Gate A discreteness={'PASS' if gate_a_pass else 'FAIL'}")
        if gate_b_pass is not None:
            cond_bits.append(f"Gate B replication={'PASS' if gate_b_pass else 'FAIL'}")
        cond = "; ".join(cond_bits) if cond_bits else "not conditioned"

        if under:
            interp = (
                f"n={n} < {self.min_n}: the 3-way split-conformal partition is under-powered; "
                f"assignment sharpness is reported for reference only, not as a gate. "
                f"Achieved coverage {achieved:.3f} (95% CI {cov_ci})."
            )
        else:
            base_txt = (
                f"At the single 0.90 target, achieved coverage {achieved:.3f} "
                f"(95% CI {cov_ci}); assignment sharpness (singleton fraction) "
                f"{sharp:.3f} (95% CI {sharp_ci}). "
            )
            gate_txt = (
                "This quantifies how sharply patients are assigned GIVEN a partition; "
                "it is meaningful only when that partition is itself certified. "
                f"Conditioning: {cond}."
            )
            interp = base_txt + gate_txt

        return GateCReframedResult(
            tumor=tumor,
            n=int(n),
            k=int(k),
            conditional_on=cond,
            target_coverage=self.primary_coverage,
            achieved_coverage=round(achieved, 4),
            coverage_ci95=cov_ci,
            coverage_achieved=cov_ok,
            sharpness_confident_fraction=round(sharp, 4),
            sharpness_ci95=sharp_ci,
            mean_set_size=res.mean_set_size[key],
            n_splits_used=int(len(cov_splits)),
            under_powered=under,
            interpretation=interp,
        )
