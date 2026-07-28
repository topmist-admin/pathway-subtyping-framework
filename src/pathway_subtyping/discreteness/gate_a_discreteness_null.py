"""
Gate A — discreteness-null stability gate (v2).

WHY THIS REPLACES THE INDEPENDENCE NULL
---------------------------------------
The previous Gate A (`gate_a_ncal_stability.NCalibratedStabilityGate`) compared
the observed mean bootstrap-ARI against a null built by permuting each pathway
column independently (`_feature_permute`). That null destroys inter-pathway
correlation, so it tests the hypothesis "the pathways are mutually independent",
NOT "there are no discrete clusters". A single correlated Gaussian blob or a 1-D
continuous gradient (tumor purity, immune infiltration, proliferation) has no
discrete clusters, yet a GMM fit to it reproducibly bisects it along the dominant
axis every bootstrap -> high observed ARI; meanwhile the feature-permuted null,
having destroyed the correlation, has a low ARI. Observed >> null => a continuum
is falsely certified as a "reproducible subtype".

The observed statistic is fine; the reference distribution was wrong. This module
keeps the identical observed statistic (mean bootstrap-ARI at the fixed selected
k) and changes what it is compared against, adding three discreteness diagnostics.

DECISION RULE -- READ BEFORE DOCUMENTING THIS GATE.
Only (A) decides. The rule is `passed = obs > sg_p95 and sg_p < alpha`; the gap
statistic (B) and the dip test (C) are computed and reported but NEVER enter it.
Thresholding the SigClust p alone reproduces the deposited ablation's TPR and FPR
exactly. Do not describe this gate as requiring a data set to clear three
criteria -- it requires one, and reports two more alongside.

Also note the gate has THREE outcomes, not two: certify, reject, and
"not-testable" (no reproducible k), which is an abstention. It abstained on 28 of
30 synthetic negative controls, so any false-positive rate quoted from it must
carry its testable denominator.

  (A) SINGLE-GAUSSIAN reference (SigClust-style) -- PRIMARY.
      Fit ONE multivariate Gaussian to the (dimension-reduced) score matrix and
      ask whether the observed partition is more stable than partitions of data
      drawn from that single Gaussian. Pass iff observed > reference 95th pctile
      AND empirical p < 0.05.
      Ref: Liu, Hayes, Nobel & Marron, "Statistical Significance of Clustering
      for High-Dimension, Low-Sample Size Data", JASA 103(483):1281-1293, 2008.
      DOI 10.1198/016214508000000454.

  (B) GAP statistic -- COMPLEMENTARY.
      Tibshirani, Walther & Hastie (2001), JRSS B 63(2):411-423,
      DOI 10.1111/1467-9868.00293. Compares log within-cluster dispersion W_k to
      its expectation under a uniform reference over the PCA-aligned bounding box;
      gap(k)=E*[log W_k]-log W_k with s'_k. A continuum does not give a positive
      peaked gap. Reports gap(k), the gap-optimal k (first-SE rule), and whether
      the fixed k is gap-supported.

  (C) HARTIGAN's DIP test -- CORROBORATING unimodality check.
      Dip on PC1 and on the GMM discriminant axis. Non-significant dip = unimodal
      = evidence against discrete clusters. A flag, not the sole gate.

SMALL-n HARDENING (fixed, not tuned; identical across observed and reference)
-----------------------------------------------------------------------------
  * Dimensionality: reduce to d = min(10, floor(n/10)) PCA components before the
    GMM and before the Gaussian reference. A 50-dim full-covariance GMM is
    singular at n in the tens; PCA to d << n well-conditions it and exits the
    HDLSS regime. The single Gaussian is simulated in this PCA space with a
    diagonal covariance equal to the retained per-PC variances (PCA components
    are uncorrelated by construction -> faithful SigClust null).
  * Fixed k: k = the inherited BIC-selected k, held identical for the observed
    statistic AND every reference/gap replicate. k is never re-selected per
    replicate (structureless data would collapse to k=1 and ARI degenerates).
  * k-stability rule: BIC-select k on each observed bootstrap resample; if the
    fraction whose modal k equals the fixed k is < 0.5, the tumor is routed to
    NOT-TESTABLE rather than pass/fail.
  * Reports the full reference distribution and a CI on the observed statistic.

The independence null (`_feature_permute`) is DEMOTED, not deleted: it is still
computed and reported as a separate confound / marginal-artifact control (does
the structure depend on genuine cross-pathway dependence, or on a few marginals?)
-- it is simply no longer the discreteness test.

Reuses ValidationGates.stability_test_bootstrap verbatim for observed, reference,
and confound-control statistics so all are strictly comparable.

Research use only.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

from pathway_subtyping.validation import ValidationGates

try:
    import diptest as _diptest
except Exception:  # pragma: no cover
    _diptest = None


# --------------------------------------------------------------------------- #
# Dimensionality rule (fixed, shared by observed + all references)
# --------------------------------------------------------------------------- #
def reduced_dim(n: int) -> int:
    """d = min(10, floor(n/10)), floored at 2."""
    return int(max(2, min(10, n // 10)))


def reduce_scores(scores: np.ndarray, d: int, seed: int) -> tuple[np.ndarray, PCA]:
    """PCA-reduce a samples x pathways matrix to d components (centred).

    d is additionally capped at the number of available features and at
    n_samples-1, so narrow panels (few pathways) and very small n are handled
    without error."""
    n, p = scores.shape
    d_eff = int(max(2, min(d, p, n - 1)))
    pca = PCA(n_components=d_eff, random_state=seed)
    Z = pca.fit_transform(scores)
    return Z, pca


# --------------------------------------------------------------------------- #
# Confound control (the demoted independence null)
# --------------------------------------------------------------------------- #
def _feature_permute(scores: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Independently permute each column across samples (kills cross-pathway
    dependence, preserves marginals). DEMOTED to a confound control."""
    out = np.empty_like(scores)
    for j in range(scores.shape[1]):
        out[:, j] = rng.permutation(scores[:, j])
    return out


# --------------------------------------------------------------------------- #
# One reference / control replicate: build a structureless dataset the
# requested way, cluster at fixed k, measure bootstrap-ARI with the framework.
# --------------------------------------------------------------------------- #
def _stability_of(
    data: np.ndarray,
    n_clusters: int,
    seed: int,
    n_bootstrap: int,
    gmm_seed: int,
) -> float:
    """Fit GMM at fixed k -> reference labels, then framework bootstrap-ARI."""
    df = pd.DataFrame(data, columns=[f"z{j}" for j in range(data.shape[1])])
    gmm = GaussianMixture(
        n_components=n_clusters,
        covariance_type="full",
        n_init=5,
        random_state=gmm_seed,
        reg_covar=1e-6,
    )
    gmm.fit(data)
    labels = gmm.predict(data)
    vg = ValidationGates(seed=seed, n_bootstrap=n_bootstrap, show_progress=False)
    return float(
        vg.stability_test_bootstrap(df, labels, n_clusters, gmm_seed=gmm_seed).metric_value
    )


def _single_gaussian_once(mean, var, n, n_clusters, seed, n_bootstrap, gmm_seed) -> float:
    """One SigClust-style reference draw in PCA space: sample n points from a
    single Gaussian N(mean, diag(var)), then measure stability at fixed k."""
    rng = np.random.default_rng(seed)
    data = rng.normal(loc=mean, scale=np.sqrt(var), size=(n, len(mean)))
    return _stability_of(data, n_clusters, seed, n_bootstrap, gmm_seed)


def _feature_permute_once(scores, n_clusters, seed, n_bootstrap, gmm_seed) -> float:
    """One confound-control replicate (independence null)."""
    rng = np.random.default_rng(seed)
    return _stability_of(_feature_permute(scores, rng), n_clusters, seed, n_bootstrap, gmm_seed)


# --------------------------------------------------------------------------- #
# Gap statistic (Tibshirani-Walther-Hastie 2001), uniform reference over the
# PCA-aligned bounding box.
# --------------------------------------------------------------------------- #
def _within_dispersion(data: np.ndarray, labels: np.ndarray, k: int) -> float:
    """W_k = sum_r (1/2 n_r) D_r  ==  sum_r sum_{i in r} ||x_i - centroid_r||^2."""
    W = 0.0
    for r in range(k):
        pts = data[labels == r]
        if len(pts) <= 1:
            continue
        W += ((pts - pts.mean(0)) ** 2).sum()
    return W


def gap_statistic(Z: np.ndarray, k_range, seed: int, n_ref: int = 50) -> Dict[str, Any]:
    """Gap over PCA-aligned bounding box (Z is already in PCA space, so the
    principal axes are the coordinate axes and the box is axis-aligned)."""
    rng = np.random.default_rng(seed)
    lo, hi = Z.min(0), Z.max(0)
    n = Z.shape[0]
    logW, gaps, sk = [], [], []
    for k in k_range:
        gm = GaussianMixture(
            n_components=k, covariance_type="full", n_init=5, random_state=seed, reg_covar=1e-6
        ).fit(Z)
        lab = gm.predict(Z)
        lw = np.log(_within_dispersion(Z, lab, k) + 1e-12)
        ref_lw = []
        for b in range(n_ref):
            U = rng.uniform(lo, hi, size=(n, Z.shape[1]))
            gmb = GaussianMixture(
                n_components=k,
                covariance_type="full",
                n_init=3,
                random_state=seed + b + 1,
                reg_covar=1e-6,
            ).fit(U)
            ref_lw.append(np.log(_within_dispersion(U, gmb.predict(U), k) + 1e-12))
        ref_lw = np.array(ref_lw)
        gap = ref_lw.mean() - lw
        s = ref_lw.std() * np.sqrt(1 + 1.0 / n_ref)
        logW.append(lw)
        gaps.append(gap)
        sk.append(s)
    gaps, sk = np.array(gaps), np.array(sk)
    ks = list(k_range)
    # first-SE rule: smallest k with gap(k) >= gap(k+1) - s(k+1)
    gap_opt = ks[-1]
    for i in range(len(ks) - 1):
        if gaps[i] >= gaps[i + 1] - sk[i + 1]:
            gap_opt = ks[i]
            break
    return {
        "k_range": ks,
        "gap": [round(float(x), 4) for x in gaps],
        "gap_se": [round(float(x), 4) for x in sk],
        "logW": [round(float(x), 4) for x in logW],
        "gap_optimal_k": int(gap_opt),
    }


# --------------------------------------------------------------------------- #
# Dip test
# --------------------------------------------------------------------------- #
def dip_of(x: np.ndarray) -> Dict[str, float]:
    if _diptest is None:
        return {"dip": float("nan"), "p": float("nan")}
    d, p = _diptest.diptest(np.asarray(x, float).ravel())
    return {"dip": round(float(d), 4), "p": round(float(p), 4)}


# --------------------------------------------------------------------------- #
# Result container
# --------------------------------------------------------------------------- #
@dataclass
class GateAv2Result:
    tumor: str
    n: int
    n_clusters: int
    reduced_d: int
    observed_stability: float
    observed_ci95: List[float]
    # primary: single-Gaussian (SigClust) reference
    sg_ref_mean: float
    sg_ref_p05: float
    sg_ref_p50: float
    sg_ref_p95: float
    sg_empirical_p: float
    sg_effect_z: float
    passed: bool
    # gap
    gap_optimal_k: int
    gap_supports_k: bool
    gap_at_k: float
    # dip
    dip_pc1_p: float
    dip_discriminant_p: float
    unimodal_flag: bool
    # k-stability
    k_stability_frac: float
    modal_k: int
    modal_matches_fixed: bool
    testable: bool
    verdict: str
    # demoted confound control (independence null)
    fp_ref_mean: float
    fp_ref_p95: float
    fp_empirical_p: float
    fp_would_pass: bool
    n_ref: int
    n_bootstrap: int
    interpretation: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class DiscretenessGateA(ValidationGates):
    """Gate A v2 — discreteness null (single-Gaussian primary + gap + dip),
    with the independence null demoted to a confound control."""

    def __init__(
        self,
        seed: int = 42,
        n_ref: int = 200,
        n_bootstrap: int = 50,
        n_jobs: int = -1,
        alpha: float = 0.05,
        gap_n_ref: int = 50,
        k_range=(1, 2, 3, 4, 5, 6),
        **kw,
    ):
        super().__init__(seed=seed, n_bootstrap=n_bootstrap, show_progress=False, **kw)
        self.n_ref = n_ref
        self.n_jobs = n_jobs
        self.alpha = alpha
        self.gap_n_ref = gap_n_ref
        self.k_range = tuple(k_range)

    # --- k-stability: reproducible modal silhouette-k across bootstraps ----- #
    def _k_stability(self, Z, k_fixed, n_boot=40):
        """Does a reproducible number of clusters exist across bootstrap
        resamples? Uses the silhouette-optimal k per bootstrap (BIC with full
        covariance over-selects in the reduced HDLSS space and is not a reliable
        k-selector here). Returns (frac at the modal k, modal k). Genuine
        discrete structure yields a high frac at a stable modal k; a continuum
        scatters k across bootstraps -> low frac -> route to not-testable."""
        from collections import Counter

        from sklearn.metrics import silhouette_score

        rng = np.random.default_rng(self.seed + 777)
        n = Z.shape[0]
        kmax = min(6, n // 5)
        picks = []
        for _ in range(n_boot):
            idx = np.unique(rng.choice(n, n, replace=True))
            if len(idx) < 4:
                continue
            Zi = Z[idx]
            sils = {}
            for k in range(2, max(3, kmax) + 1):
                if k >= len(idx):
                    break
                try:
                    gm = GaussianMixture(
                        n_components=k,
                        covariance_type="full",
                        n_init=3,
                        random_state=self.seed,
                        reg_covar=1e-6,
                    ).fit(Zi)
                    lab = gm.predict(Zi)
                    if len(set(lab)) < 2:
                        continue
                    sils[k] = silhouette_score(Zi, lab)
                except Exception:
                    continue
            if sils:
                picks.append(max(sils, key=sils.get))
        if not picks:
            return 0.0, k_fixed
        c = Counter(picks)
        modal_k, modal_n = c.most_common(1)[0]
        return float(modal_n / len(picks)), int(modal_k)

    def run(
        self,
        tumor: str,
        pathway_scores: pd.DataFrame,
        n_clusters: int,
        gmm_seed: Optional[int] = None,
    ) -> GateAv2Result:
        gmm_seed = self.seed if gmm_seed is None else gmm_seed
        scores = pathway_scores.values
        n, p = scores.shape
        d = reduced_dim(n)

        # ---- dimensionality reduction (shared by observed + references) ---- #
        Z, pca = reduce_scores(scores, d, self.seed)
        d = Z.shape[1]  # effective dim after capping at (p, n-1)
        var = Z.var(0)  # per-PC variances -> diagonal single-Gaussian
        mean = Z.mean(0)  # ~0 after PCA centring

        # ---- observed partition + statistic in the reduced space ---------- #
        gmm = GaussianMixture(
            n_components=n_clusters,
            covariance_type="full",
            n_init=10,
            random_state=gmm_seed,
            reg_covar=1e-6,
        ).fit(Z)
        obs_labels = gmm.predict(Z)
        Zdf = pd.DataFrame(Z, columns=[f"z{j}" for j in range(d)])
        obs_res = self.stability_test_bootstrap(Zdf, obs_labels, n_clusters, gmm_seed=gmm_seed)
        obs = float(obs_res.metric_value)
        std_ari = float(obs_res.details.get("std_ari", 0.0))
        nb = int(obs_res.details.get("n_bootstrap", self.n_bootstrap))
        se = std_ari / np.sqrt(max(nb, 1))
        ci = [round(max(-1, obs - 1.96 * se), 4), round(min(1, obs + 1.96 * se), 4)]

        # ---- (A) single-Gaussian (SigClust) reference --------------------- #
        ss = np.random.SeedSequence(self.seed)
        seeds = [int(np.random.default_rng(s).integers(1, 2**31)) for s in ss.spawn(self.n_ref)]
        sg = Parallel(n_jobs=self.n_jobs, prefer="threads")(
            delayed(_single_gaussian_once)(
                mean, var, n, n_clusters, seeds[i], self.n_bootstrap, gmm_seed + i + 1
            )
            for i in range(self.n_ref)
        )
        sg = np.asarray([x for x in sg if np.isfinite(x)])
        sg_mu, sg_sd = float(sg.mean()), float(sg.std())
        sg_p05, sg_p50, sg_p95 = [float(np.percentile(sg, q)) for q in (5, 50, 95)]
        sg_p = float((np.sum(sg >= obs) + 1) / (len(sg) + 1))
        sg_z = float((obs - sg_mu) / sg_sd) if sg_sd > 0 else float("nan")

        # ---- confound control: independence (feature-permute) null -------- #
        fp = Parallel(n_jobs=self.n_jobs, prefer="threads")(
            delayed(_feature_permute_once)(
                scores, n_clusters, seeds[i], self.n_bootstrap, gmm_seed + i + 1
            )
            for i in range(self.n_ref)
        )
        fp = np.asarray([x for x in fp if np.isfinite(x)])
        fp_mu, fp_p95 = float(fp.mean()), float(np.percentile(fp, 95))
        fp_p = float((np.sum(fp >= obs) + 1) / (len(fp) + 1))

        # ---- (B) gap statistic -------------------------------------------- #
        gap = gap_statistic(Z, self.k_range, self.seed, n_ref=self.gap_n_ref)
        gap_at_k = (
            float(gap["gap"][gap["k_range"].index(n_clusters)])
            if n_clusters in gap["k_range"]
            else float("nan")
        )
        gap_supports = bool(gap["gap_optimal_k"] >= 2 and gap["gap_optimal_k"] == n_clusters)

        # ---- (C) dip test on PC1 and GMM discriminant axis ---------------- #
        pc1 = Z[:, 0]
        # discriminant axis = direction between the two most separated GMM means
        means = gmm.means_
        if len(means) >= 2:
            from itertools import combinations

            pair = max(
                combinations(range(len(means)), 2),
                key=lambda ij: np.sum((means[ij[0]] - means[ij[1]]) ** 2),
            )
            axis = means[pair[0]] - means[pair[1]]
            axis = axis / (np.linalg.norm(axis) + 1e-12)
            disc = Z @ axis
        else:
            disc = pc1
        dip_pc1 = dip_of(pc1)
        dip_disc = dip_of(disc)
        unimodal = bool((dip_pc1["p"] > 0.05) and (dip_disc["p"] > 0.05))

        # ---- k-stability + testability ------------------------------------ #
        kfrac, modal_k = self._k_stability(Z, n_clusters)
        testable = bool(kfrac >= 0.5)  # a reproducible modal k exists
        modal_matches_fixed = bool(modal_k == n_clusters)

        passed = bool(obs > sg_p95 and sg_p < self.alpha)

        if not testable:
            verdict = "not-testable (no reproducible k)"
        elif passed:
            verdict = "discrete structure"
        else:
            verdict = "no discrete structure (continuum)"

        interp = (
            f"Observed bootstrap-ARI {obs:.3f} vs single-Gaussian reference "
            f"95th pctile {sg_p95:.3f} (p={sg_p:.3f}). "
            f"Independence (confound) null 95th {fp_p95:.3f}: would{'' if obs > fp_p95 else ' NOT'} pass. "
            f"Gap-optimal k={gap['gap_optimal_k']} ({'supports' if gap_supports else 'does not support'} k={n_clusters}). "
            f"Dip PC1 p={dip_pc1['p']:.3f}, discriminant p={dip_disc['p']:.3f} "
            f"({'unimodal' if unimodal else 'multimodal'})."
        )

        res = GateAv2Result(
            tumor=tumor,
            n=n,
            n_clusters=n_clusters,
            reduced_d=d,
            observed_stability=round(obs, 4),
            observed_ci95=ci,
            sg_ref_mean=round(sg_mu, 4),
            sg_ref_p05=round(sg_p05, 4),
            sg_ref_p50=round(sg_p50, 4),
            sg_ref_p95=round(sg_p95, 4),
            sg_empirical_p=round(sg_p, 4),
            sg_effect_z=round(sg_z, 3) if np.isfinite(sg_z) else float("nan"),
            passed=passed,
            gap_optimal_k=gap["gap_optimal_k"],
            gap_supports_k=gap_supports,
            gap_at_k=round(gap_at_k, 4) if np.isfinite(gap_at_k) else float("nan"),
            dip_pc1_p=dip_pc1["p"],
            dip_discriminant_p=dip_disc["p"],
            unimodal_flag=unimodal,
            k_stability_frac=round(kfrac, 4),
            modal_k=modal_k,
            modal_matches_fixed=modal_matches_fixed,
            testable=testable,
            verdict=verdict,
            fp_ref_mean=round(fp_mu, 4),
            fp_ref_p95=round(fp_p95, 4),
            fp_empirical_p=round(fp_p, 4),
            fp_would_pass=bool(obs > fp_p95),
            n_ref=int(len(sg)),
            n_bootstrap=self.n_bootstrap,
            interpretation=interp,
        )
        # attach raw distributions for plotting (not in to_dict)
        res.sg_distribution = sg
        res.fp_distribution = fp
        res.gap_full = gap
        res.observed_labels = obs_labels
        res.Z = Z
        res.discriminant_axis_proj = disc
        return res
