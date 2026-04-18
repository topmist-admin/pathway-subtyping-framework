"""
Invariant Causal Prediction (ICP) for pathway-level relationships.

ICP identifies the subset of candidate-parent variables whose regression
on the target produces residuals that are *invariant* across
environments — the hallmark of causal direction. The identifiable
causal parents are the intersection of all invariant subsets.

Reference:
    Peters, Buhlmann, Meinshausen (2016). "Causal inference using
    invariant prediction: identification and confidence intervals."
    JRSS-B 78(5).

This implementation is intentionally light-weight:

    - Linear regression per subset, test whether residual means match
      across environments via a one-way ANOVA F-test.
    - Small candidate sets only (up to ~6 candidate parents); larger
      problems should use `dowhy` or a dedicated ICP package.
    - Deterministic given a fixed random seed.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import itertools
import logging
import math
from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import f as f_distribution

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Core statistics
# --------------------------------------------------------------------------- #

def invariance_pvalue(
    residuals: np.ndarray,
    environments: np.ndarray,
) -> float:
    """Combined mean+variance invariance p-value across environments.

    Returns the minimum of two tests:

        * One-way ANOVA on residual **means** across environments.
        * Levene-style test on residual **variances** across environments.

    This combined criterion is robust against the common ICP failure
    mode where a spurious subset (e.g., including a *child* of the
    target) drives residual means to zero in every environment but
    leaves environment-specific variance intact. Either rejection
    shrinks the combined p-value and prevents such a subset from being
    declared invariant.

    Under the null (residuals share mean *and* variance across envs),
    both tests concentrate near 1; under the alternative, at least one
    test gets small.

    Args:
        residuals: 1D residual vector.
        environments: Integer / string environment label per observation.

    Returns:
        p-value in [0, 1]. 1.0 when a single environment is supplied or
        residuals are trivially constant.
    """
    residuals = np.asarray(residuals, dtype=float)
    environments = np.asarray(environments)

    unique_envs = np.unique(environments)
    k = len(unique_envs)
    n = len(residuals)
    if k < 2 or n <= k:
        return 1.0

    # --- Mean test (ANOVA) ----------------------------------------------
    grand_mean = residuals.mean()
    ss_between = 0.0
    ss_within = 0.0
    abs_dev_groups: List[np.ndarray] = []
    for env in unique_envs:
        mask = environments == env
        block = residuals[mask]
        if len(block) == 0:
            continue
        env_mean = block.mean()
        ss_between += len(block) * (env_mean - grand_mean) ** 2
        ss_within += np.sum((block - env_mean) ** 2)
        abs_dev_groups.append(np.abs(block - env_mean))

    df_between = k - 1
    df_within = n - k
    if df_within <= 0 or ss_within <= 0.0:
        return 1.0
    ms_within_mean = ss_within / df_within
    if ms_within_mean == 0:
        p_mean = 1.0
    else:
        f_mean = (ss_between / df_between) / ms_within_mean
        p_mean = float(f_distribution.sf(f_mean, df_between, df_within))

    # --- Variance test (Levene, on abs deviations from per-env median)---
    # Equivalent to ANOVA on |residual - env_mean|.
    abs_dev = np.concatenate(abs_dev_groups)
    grand_abs_mean = abs_dev.mean()
    ss_between_v = 0.0
    ss_within_v = 0.0
    offset = 0
    for block_abs in abs_dev_groups:
        m = block_abs.mean()
        ss_between_v += len(block_abs) * (m - grand_abs_mean) ** 2
        ss_within_v += np.sum((block_abs - m) ** 2)
    if ss_within_v <= 0.0:
        p_var = 1.0
    else:
        ms_within_var = ss_within_v / df_within
        if ms_within_var == 0:
            p_var = 1.0
        else:
            f_var = (ss_between_v / df_between) / ms_within_var
            p_var = float(f_distribution.sf(f_var, df_between, df_within))

    return float(min(p_mean, p_var))


def _fit_residuals(
    X: np.ndarray,
    y: np.ndarray,
    ridge: float = 1e-6,
) -> np.ndarray:
    """Linear fit residuals with an intercept column; tiny ridge for stability."""
    if X.size == 0:
        return y - y.mean()
    Xa = np.hstack([X, np.ones((len(X), 1))])
    d = Xa.shape[1]
    reg = ridge * np.eye(d)
    reg[-1, -1] = 0.0
    beta = np.linalg.solve(Xa.T @ Xa + reg, Xa.T @ y)
    return y - Xa @ beta


# --------------------------------------------------------------------------- #
# Result types
# --------------------------------------------------------------------------- #

@dataclass
class CausalParentReport:
    """Outcome of an :class:`InvariantPathwayPredictor` run.

    Attributes:
        target: The target variable.
        candidate_parents: Candidate parent names searched over.
        invariant_subsets: Subsets whose residuals passed the invariance
            test at ``alpha``.
        identifiable_parents: Intersection of ``invariant_subsets`` —
            the features that appear in *every* invariant subset, which
            under ICP assumptions are the identifiable causal parents.
        alpha: Significance level used for the invariance test.
        subset_pvalues: Dict mapping each tested subset to its ANOVA
            p-value on residual invariance.
    """

    target: str
    candidate_parents: List[str]
    invariant_subsets: List[FrozenSet[str]]
    identifiable_parents: FrozenSet[str]
    alpha: float
    subset_pvalues: Dict[FrozenSet[str], float] = field(default_factory=dict)

    @property
    def n_invariant_subsets(self) -> int:
        return len(self.invariant_subsets)

    def recall_against(self, ground_truth: Iterable[str]) -> float:
        """Fraction of ``ground_truth`` parents identified."""
        gt = set(ground_truth)
        if not gt:
            return float("nan")
        return float(len(self.identifiable_parents & gt)) / float(len(gt))

    def precision_against(self, ground_truth: Iterable[str]) -> float:
        """Fraction of identified parents that are in ``ground_truth``."""
        if not self.identifiable_parents:
            return float("nan")
        gt = set(ground_truth)
        return float(len(self.identifiable_parents & gt)) / float(
            len(self.identifiable_parents)
        )

    def summary(self) -> str:
        return (
            f"CausalParentReport(target={self.target!r}, "
            f"identifiable={sorted(self.identifiable_parents)}, "
            f"invariant_subsets={self.n_invariant_subsets}, "
            f"alpha={self.alpha})"
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "target": self.target,
            "candidate_parents": list(self.candidate_parents),
            "identifiable_parents": sorted(self.identifiable_parents),
            "n_invariant_subsets": self.n_invariant_subsets,
            "invariant_subsets": [sorted(s) for s in self.invariant_subsets],
            "alpha": self.alpha,
        }


# --------------------------------------------------------------------------- #
# ICP driver
# --------------------------------------------------------------------------- #

class InvariantPathwayPredictor:
    """Fit ICP for a single target variable over pathway-level data.

    Usage::

        predictor = InvariantPathwayPredictor(alpha=0.05, max_subset_size=3)
        report = predictor.fit(
            X=pathway_scores,      # samples x pathways
            y=target_series,       # one of the pathway columns
            target_name="PATH_0",
            environments=env_labels,
        )
        print(report.identifiable_parents)
    """

    def __init__(
        self,
        alpha: float = 0.05,
        max_subset_size: Optional[int] = 4,
        include_empty_set: bool = True,
    ) -> None:
        if not 0.0 < alpha < 1.0:
            raise ValueError("alpha must be in (0, 1)")
        self.alpha = float(alpha)
        self.max_subset_size = max_subset_size
        self.include_empty_set = bool(include_empty_set)

    # -------------------------------------------------------------- fit ---
    def fit(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        target_name: str,
        environments: Sequence[Any],
    ) -> CausalParentReport:
        """Search invariant subsets and return the identifiable parents.

        Args:
            X: Samples × candidate-parents DataFrame. Must NOT contain
                the target column (caller's responsibility).
            y: Target values (same index as ``X``).
            target_name: Name of the target variable, for reporting.
            environments: Environment label per sample; any hashable
                value is accepted.
        """
        if len(X) != len(y) or len(X) != len(environments):
            raise ValueError("X, y, environments must share the same length")
        if target_name in X.columns:
            raise ValueError(
                f"target_name={target_name!r} must not be a column of X"
            )

        X_values = X.to_numpy(dtype=float)
        y_values = np.asarray(y, dtype=float)
        env_values = np.asarray(environments)

        candidate_parents = list(X.columns)
        if self.max_subset_size is None:
            max_k = len(candidate_parents)
        else:
            max_k = min(int(self.max_subset_size), len(candidate_parents))

        # Enumerate subsets
        subsets: List[FrozenSet[str]] = []
        if self.include_empty_set:
            subsets.append(frozenset())
        for k in range(1, max_k + 1):
            for combo in itertools.combinations(candidate_parents, k):
                subsets.append(frozenset(combo))

        invariant: List[FrozenSet[str]] = []
        pvals: Dict[FrozenSet[str], float] = {}
        col_index = {c: i for i, c in enumerate(candidate_parents)}

        for subset in subsets:
            if not subset:
                X_sub = np.empty((len(X_values), 0))
            else:
                cols = [col_index[p] for p in subset]
                X_sub = X_values[:, cols]
            residuals = _fit_residuals(X_sub, y_values)
            p = invariance_pvalue(residuals, env_values)
            pvals[subset] = p
            if p > self.alpha:
                invariant.append(subset)

        # Identifiable parents: intersection of all invariant subsets
        if not invariant:
            identifiable: FrozenSet[str] = frozenset()
        else:
            identifiable = frozenset.intersection(*invariant)

        logger.info(
            "[ICP] target=%s invariant_subsets=%d identifiable=%s",
            target_name, len(invariant), sorted(identifiable),
        )

        return CausalParentReport(
            target=target_name,
            candidate_parents=candidate_parents,
            invariant_subsets=invariant,
            identifiable_parents=identifiable,
            alpha=self.alpha,
            subset_pvalues=pvals,
        )
