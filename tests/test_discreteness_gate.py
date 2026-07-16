"""
Synthetic-control validation for the discreteness-aware Gate A (v0.8.0).

These are the controls that demonstrate the gate does what it claims: it must
FAIL structureless data and continuous gradients (which the old independence
null wrongly passed) and PASS genuinely separated clusters. Ground truth is
known by construction, so the asserts are on the verdict direction, not on any
real-data result.

Reduced reference/bootstrap counts keep the run fast while preserving the
direction; seeds are fixed for determinism.
"""
import numpy as np
import pandas as pd
import pytest

from pathway_subtyping.discreteness import DiscretenessGateA, ReframedMembershipGate


SEED = 42
N = 48          # small-n regime (tens of samples), like a rare-tumor cohort
P = 8           # pathway-score dimensions


def _df(mat: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(mat, columns=[f"pw{j}" for j in range(mat.shape[1])])


def _single_gaussian(n=N, p=P, seed=SEED):
    """No clusters: one isotropic Gaussian blob."""
    rng = np.random.default_rng(seed)
    return _df(rng.normal(size=(n, p)))


def _continuum(n=N, p=P, seed=SEED):
    """No clusters: a 1-D gradient (points along a line + small orthogonal noise).
    This is the case the OLD independence null falsely certified as a subtype."""
    rng = np.random.default_rng(seed)
    t = np.linspace(-3, 3, n)
    direction = rng.normal(size=p)
    direction /= np.linalg.norm(direction)
    mat = np.outer(t, direction) + 0.35 * rng.normal(size=(n, p))
    return _df(mat)


def _two_clusters(n=N, p=P, seed=SEED, sep=5.0):
    """Genuine discrete structure: two well-separated Gaussian clusters."""
    rng = np.random.default_rng(seed)
    half = n // 2
    offset = np.zeros(p)
    offset[0] = sep
    a = rng.normal(size=(half, p)) - offset / 2
    b = rng.normal(size=(n - half, p)) + offset / 2
    return _df(np.vstack([a, b]))


def _gate():
    # reduced counts for test speed; direction of the verdict is robust to N_ref
    return DiscretenessGateA(seed=SEED, n_ref=40, n_bootstrap=12,
                             n_jobs=1, gap_n_ref=12, k_range=(1, 2, 3))


@pytest.mark.slow
def test_discreteness_fails_single_gaussian():
    res = _gate().run("neg_gaussian", _single_gaussian(), n_clusters=2)
    assert res.passed is False
    assert "discrete structure" not in res.verdict or res.verdict.startswith("no")


@pytest.mark.slow
def test_discreteness_fails_continuum():
    # the decisive case: the OLD independence null passed this; the NEW must fail it
    res = _gate().run("continuum", _continuum(), n_clusters=2)
    assert res.passed is False


@pytest.mark.slow
def test_discreteness_passes_two_clusters():
    res = _gate().run("two_clusters", _two_clusters(), n_clusters=2)
    assert res.passed is True
    assert res.verdict == "discrete structure"
    # sanity: observed stability clears the single-Gaussian reference 95th pctile
    assert res.observed_stability > res.sg_ref_p95


@pytest.mark.slow
def test_reframed_gate_c_runs_and_is_conditional():
    """Smoke test: reframed Gate C reports a single-target coverage with a CI and
    carries the explicit conditional-on-A/B framing."""
    from sklearn.mixture import GaussianMixture
    scores = _two_clusters()
    labels = GaussianMixture(n_components=2, random_state=SEED).fit_predict(scores.values)
    res = ReframedMembershipGate(seed=SEED, n_splits=20).run(
        "two_clusters", scores, labels, k=2, gate_a_pass=True, gate_b_pass=True)
    assert res.target_coverage == 0.90
    assert np.isfinite(res.achieved_coverage)
    assert len(res.coverage_ci95) == 2
    assert "Gate A discreteness=PASS" in res.conditional_on
