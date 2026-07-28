"""
Uncertainty Quantification Layer for the Pathway Subtyping Framework.

Provides three complementary uncertainty channels plus a calibration
assessment toolkit:

Submodules:
    - conformal: Split-conformal prediction intervals (distribution-free)
    - bootstrap: Non-parametric bootstrap over cells for MSV scores
    - bayesian_pathway: Bayesian variant of pathway_gmm with posterior sampling
    - calibration: Reliability diagrams, ECE, Brier score

Every MSV output today (maturation_score, stress_score, drift_score, etc.)
returns a single scalar. A researcher reading stress_score = 0.72 cannot
distinguish "confidently in a stress state" from "model is uncertain and
true value could be anywhere from 0.4 to 0.9." This module closes that gap.

Research use only. Not for clinical decision-making.
"""

from .bayesian_pathway import BayesianPathwayGMM, PosteriorSample
from .bootstrap import BootstrapMSV, BootstrapResult
from .calibration import CalibrationReport, brier_score, ece, reliability_curve
from .conformal import ConformalInterval, ConformalPathwayPredictor

__all__ = [
    "BootstrapMSV",
    "BootstrapResult",
    "CalibrationReport",
    "ConformalInterval",
    "ConformalPathwayPredictor",
    "BayesianPathwayGMM",
    "PosteriorSample",
    "brier_score",
    "ece",
    "reliability_curve",
]
