"""
Gate C — Conformal (calibrated-uncertainty) subtype-membership confidence.

Extends the pathway-subtyping framework's split-conformal predictor
(ConformalPathwayPredictor, classification/APS mode) and the existing
conformal-coverage routine to emit a PER-PATIENT calibrated confidence on
subtype membership.

Method
------
The subtype model is a Gaussian mixture fit on the cohort's pathway scores; its
posterior P(subtype | patient) is the classification score function. Using a
held-out calibration split, split-conformal (APS-style nonconformity 1 - p_true)
computes a quantile that guarantees marginal coverage 1 - alpha. Each patient
then receives a conformal PREDICTION SET of candidate subtypes:

    * singleton set  → confidently assigned to one subtype
    * multi-label set → ambiguous (cannot be confidently assigned)
    * (empty sets are promoted to the top-1 label)

Reported per cohort, at target coverages {0.80, 0.90, 0.95}, averaged over
repeated calibration/test splits (seeded):
    * empirical coverage vs target (calibration quality)
    * mean prediction-set size
    * fraction of confidently-assigned (singleton) patients
    * per-patient confident/ambiguous flag + membership probability

This mirrors the structure of the framework's existing
results/f1_validation/*_conformal_coverage.json output.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import train_test_split

from pathway_subtyping.uncertainty.conformal import ConformalPathwayPredictor


@dataclass
class GateCResult:
    tumor: str
    n: int
    k: int
    target_coverages: List[float]
    empirical_coverage: Dict[str, float]
    mean_set_size: Dict[str, float]
    confident_fraction: Dict[str, float]
    coverage_calibrated: bool
    n_splits: int
    primary_coverage: float
    primary_confident_fraction: float
    per_patient: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d.pop("per_patient", None)
        return d


class ConformalMembershipGate:
    """Gate C — per-patient conformal subtype-membership confidence."""

    def __init__(
        self,
        seed: int = 42,
        n_splits: int = 50,
        target_coverages=(0.80, 0.90, 0.95),
        cal_fraction: float = 0.5,
        primary_coverage: float = 0.90,
        calib_tol: float = 0.05,
    ):
        self.seed = seed
        self.n_splits = n_splits
        self.target_coverages = list(target_coverages)
        self.cal_fraction = cal_fraction
        self.primary_coverage = primary_coverage
        self.calib_tol = calib_tol

    def _fit_gmm(self, X, k, seed):
        gmm = GaussianMixture(
            n_components=k, covariance_type="full", n_init=10, random_state=seed, reg_covar=1e-6
        )
        gmm.fit(X)
        return gmm

    def membership_gate(
        self, tumor: str, pathway_scores: pd.DataFrame, cluster_labels: np.ndarray, k: int
    ) -> GateCResult:
        X = pathway_scores.values
        y = np.asarray(cluster_labels)
        n = X.shape[0]

        cov_acc = {f"{c:.2f}": [] for c in self.target_coverages}
        size_acc = {f"{c:.2f}": [] for c in self.target_coverages}
        conf_acc = {f"{c:.2f}": [] for c in self.target_coverages}
        # accumulate per-patient confident counts at primary coverage
        per_patient_conf = np.zeros(n)
        per_patient_seen = np.zeros(n)
        pp_key = f"{self.primary_coverage:.2f}"

        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler

        gmm_labels = np.arange(k)
        rng = np.random.default_rng(self.seed)
        for s in range(self.n_splits):
            split_seed = int(rng.integers(1, 2**31))
            idx = np.arange(n)
            # Proper 3-way split-conformal:
            #   train  -> fit the subtype CLASSIFIER (honest generalization)
            #   cal    -> nonconformity quantile
            #   test   -> evaluate coverage / set size
            # A classifier (not the generative GMM that produced the labels) is
            # used so posteriors are genuinely graded; the train/cal/test split
            # makes the calibration and evaluation out-of-training and thus honest.
            # size folds so the calibration fold has >=10 (PSF conformal floor)
            # where n allows: train ~50%, calibration >=10, test the remainder.
            n_train = max(k * 2, int(round(0.45 * n)))
            n_cal = max(10, int(round(0.30 * n)))
            if n - n_train - n_cal < 5:  # cohort too small for a clean 3-way split
                n_train = int(round(0.5 * (n - 10)))
                n_cal = 10
                if n - n_train - n_cal < 3:
                    continue
            try:
                tr_idx, rest = train_test_split(
                    idx, train_size=n_train, random_state=split_seed, stratify=y
                )
                cal_idx, test_idx = train_test_split(
                    rest, train_size=n_cal, random_state=split_seed, stratify=y[rest]
                )
            except ValueError:
                tr_idx, rest = train_test_split(idx, train_size=n_train, random_state=split_seed)
                cal_idx, test_idx = train_test_split(
                    rest, train_size=n_cal, random_state=split_seed
                )
            if len(cal_idx) < 10 or len(test_idx) < 3 or len(np.unique(y[tr_idx])) < k:
                continue

            scaler = StandardScaler().fit(X[tr_idx])
            clf = LogisticRegression(max_iter=2000, C=1.0)
            clf.fit(scaler.transform(X[tr_idx]), y[tr_idx])
            clf_classes = clf.classes_

            def score_fn(Z, _clf=clf, _sc=scaler, _cls=clf_classes):
                P = _clf.predict_proba(_sc.transform(np.asarray(Z)))
                # reindex columns to full 0..k-1 label space
                full = np.zeros((P.shape[0], k))
                for j, cl in enumerate(_cls):
                    full[:, int(cl)] = P[:, j]
                return full

            for cov in self.target_coverages:
                pred = ConformalPathwayPredictor(
                    score_fn=score_fn, coverage=cov, mode="classification"
                )
                pred.calibrate(X[cal_idx], y[cal_idx], labels=list(gmm_labels))
                intervals = pred.predict(X[test_idx])
                # empirical coverage: true label in prediction set
                hits, sizes, singles = 0, [], 0
                for iv, yi in zip(intervals, y[test_idx]):
                    ls = iv.label_set if iv.label_set else [int(np.round(iv.prediction))]
                    if int(yi) in [int(x) for x in ls]:
                        hits += 1
                    sizes.append(len(ls))
                    if len(ls) == 1:
                        singles += 1
                m = len(intervals)
                key = f"{cov:.2f}"
                cov_acc[key].append(hits / m)
                size_acc[key].append(float(np.mean(sizes)))
                conf_acc[key].append(singles / m)
                if key == pp_key:
                    for j, iv in zip(test_idx, intervals):
                        ls = iv.label_set if iv.label_set else [0]
                        per_patient_seen[j] += 1
                        if len(ls) == 1:
                            per_patient_conf[j] += 1

        emp_cov = {k_: round(float(np.mean(v)), 4) for k_, v in cov_acc.items()}
        mean_size = {k_: round(float(np.mean(v)), 4) for k_, v in size_acc.items()}
        conf_frac = {k_: round(float(np.mean(v)), 4) for k_, v in conf_acc.items()}
        # retain per-split arrays (used by the reframed Gate C for a coverage CI)
        self._last_cov_per_split = {k_: list(v) for k_, v in cov_acc.items()}
        self._last_conf_per_split = {k_: list(v) for k_, v in conf_acc.items()}
        # calibration: empirical coverage within tol of every target
        calibrated = all(
            abs(emp_cov[f"{c:.2f}"] - c) <= self.calib_tol for c in self.target_coverages
        )

        # per-patient confident probability at primary coverage
        with np.errstate(invalid="ignore"):
            pp_prob = np.where(per_patient_seen > 0, per_patient_conf / per_patient_seen, np.nan)
        per_patient = {
            "sample_id": list(pathway_scores.index),
            "subtype": [int(x) for x in y],
            "confident_prob": [round(float(x), 3) if np.isfinite(x) else None for x in pp_prob],
            "confident_flag": [bool(x >= 0.5) if np.isfinite(x) else None for x in pp_prob],
        }

        return GateCResult(
            tumor=tumor,
            n=int(n),
            k=int(k),
            target_coverages=self.target_coverages,
            empirical_coverage=emp_cov,
            mean_set_size=mean_size,
            confident_fraction=conf_frac,
            coverage_calibrated=bool(calibrated),
            n_splits=self.n_splits,
            primary_coverage=self.primary_coverage,
            primary_confident_fraction=conf_frac[pp_key],
            per_patient=per_patient,
        )
