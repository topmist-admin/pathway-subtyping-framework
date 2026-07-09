"""
Guarded clustering-agreement metrics.

Background
----------
``sklearn.metrics.adjusted_rand_score`` returns a *misleading* value on
degenerate input rather than raising:

    >>> from sklearn.metrics import adjusted_rand_score
    >>> adjusted_rand_score([], [])                 # empty
    1.0
    >>> adjusted_rand_score([0, 0, 0, 0], [0, 0, 0, 0])  # 1 true cluster
    1.0
    >>> adjusted_rand_score([0], [0])               # 1 sample
    1.0

A "perfect" ARI of 1.0 from empty / single-cluster / single-sample ground
truth is an arithmetic artifact, not measured agreement. This exact failure
mode produced 14 degenerate rows in the 47-dataset calibration benchmark
(`n_true_clusters == 1`), 13 of which carried a spurious ARI = 1.0 that drove
the retracted R^2 = 0.889 adaptive-threshold model. See
``CORRECTION_2026-07/ERRATUM_2026-07-08.md`` and ``KNOWN-ISSUES.md`` (NB17).

Note the 14th degenerate row (GEO GSE136196: 1 true cluster, multi-way
prediction) returns ARI = 0.0, *not* 1.0 — so a guard keyed on the ARI *value*
(== 1.0) misses it. The guard below keys on the *ground-truth structure*
(< 2 distinct true clusters, too few samples, or empty), which catches all 14.

Use ``safe_adjusted_rand_score`` (or ``ari_with_validity`` when you need to
tabulate a machine-readable reason) everywhere a ground truth may be degenerate
— benchmark scoring, cross-cohort transfer, ground-truth validation. Bootstrap
*stability* compares two non-trivial clusterings and is not a ground-truth
comparison; it is unaffected, but the guard is safe to use there too.

Research use only. Not for clinical decision-making.
"""

import logging
import math
from typing import Optional, Tuple

import numpy as np
from sklearn.metrics import adjusted_rand_score

logger = logging.getLogger(__name__)

# Minimum requirements for a *meaningful* adjusted Rand index.
MIN_SAMPLES = 2
MIN_TRUE_CLUSTERS = 2


def _n_distinct(labels) -> int:
    """Number of distinct labels, without assuming the input array type."""
    arr = np.asarray(list(labels))
    if arr.size == 0:
        return 0
    return int(len(np.unique(arr)))


def ari_degenerate_reason(
    labels_true,
    labels_pred,
    *,
    min_samples: int = MIN_SAMPLES,
    min_true_clusters: int = MIN_TRUE_CLUSTERS,
) -> Optional[str]:
    """
    Return a short reason string if (labels_true, labels_pred) is a degenerate
    input for the adjusted Rand index, else ``None``.

    Degenerate cases (ARI is undefined / arithmetically misleading):
      * empty input;
      * fewer than ``min_samples`` samples;
      * fewer than ``min_true_clusters`` distinct *ground-truth* clusters
        (no structure to recover — the empty-input / n_true_clusters=1 artifact).

    A mismatch in length is a *caller bug*, not degeneracy; it is reported so the
    caller can raise.
    """
    n_true = len(list(labels_true))
    n_pred = len(list(labels_pred))
    if n_true != n_pred:
        return f"length_mismatch(true={n_true},pred={n_pred})"
    if n_true == 0:
        return "empty_input"
    if n_true < min_samples:
        return f"too_few_samples(n={n_true})"
    n_true_clusters = _n_distinct(labels_true)
    if n_true_clusters < min_true_clusters:
        return f"degenerate_ground_truth(n_true_clusters={n_true_clusters})"
    return None


def ari_with_validity(
    labels_true,
    labels_pred,
    *,
    min_samples: int = MIN_SAMPLES,
    min_true_clusters: int = MIN_TRUE_CLUSTERS,
) -> Tuple[float, bool, Optional[str]]:
    """
    Adjusted Rand index with an explicit validity flag and reason.

    Returns
    -------
    (value, valid, reason)
        ``value`` is the ARI (``float``) when valid, else ``float('nan')``.
        ``valid`` is ``False`` for degenerate input.
        ``reason`` is ``None`` when valid, else a short degeneracy code.

    Raises
    ------
    ValueError
        If ``labels_true`` and ``labels_pred`` differ in length (a caller bug).
    """
    reason = ari_degenerate_reason(
        labels_true,
        labels_pred,
        min_samples=min_samples,
        min_true_clusters=min_true_clusters,
    )
    if reason is not None and reason.startswith("length_mismatch"):
        raise ValueError(
            f"adjusted_rand_score requires equal-length inputs: {reason}"
        )
    if reason is not None:
        logger.warning(
            "adjusted_rand_score skipped: %s (returning NaN, not a spurious score)",
            reason,
        )
        return math.nan, False, reason
    return float(adjusted_rand_score(labels_true, labels_pred)), True, None


def safe_adjusted_rand_score(
    labels_true,
    labels_pred,
    *,
    min_samples: int = MIN_SAMPLES,
    min_true_clusters: int = MIN_TRUE_CLUSTERS,
    default: float = math.nan,
) -> float:
    """
    Drop-in replacement for ``sklearn.metrics.adjusted_rand_score`` that returns
    ``default`` (NaN by default) on degenerate input instead of a misleading
    "perfect" score.

    Raises ``ValueError`` on length mismatch (a real caller bug).
    """
    value, valid, _ = ari_with_validity(
        labels_true,
        labels_pred,
        min_samples=min_samples,
        min_true_clusters=min_true_clusters,
    )
    return value if valid else default
