"""
Molecular QC Layer for the Pathway Subtyping Framework.

Provides pathway-level quality control for manufactured and engineered cells,
moving from "pathway activity profiling" to "pathway resolution monitoring."

Workstream A (KG-independent):
    - offtarget: Off-target pathway activation detection (F6)
    - heterogeneity: Batch heterogeneity profiling (F7)
    - dosage: Pathway dosage & stoichiometry (F8)
    - atlas: Reference atlas comparison (F12)

Workstream B (KG-dependent, requires KG builder):
    - cascade: Incomplete cascade detection (F1)
    - temporal: Temporal pathway tracking (F2)
    - tension: Active tension scoring (F3)
    - gates: Resolution gates (F4)
    - drift: Drift detection (F5)
    - crosstalk: Pathway crosstalk & interference (F9)
    - feedback: Feedback loop integrity (F10)
    - stress: Environmental stress fingerprinting (F11)

Research use only. Not for clinical decision-making.
"""

# Workstream A — KG-independent
from .atlas import (
    AtlasComparator,
    AtlasComparisonResult,
    AtlasEntry,
    CellMapping,
)
from .dosage import (
    DosageAnalyzer,
    DosageAnalysisResult,
    DosageState,
    PathwayDosageResult,
    StoichiometryResult,
)
from .heterogeneity import (
    HeterogeneityProfiler,
    HeterogeneityResult,
    SubpopulationInfo,
)
from .offtarget import (
    ActivationClass,
    OffTargetDetector,
    OffTargetResult,
    PathwayActivationResult,
)

# Workstream B — KG-dependent
from .cascade import (
    CascadeAnalyzer,
    CascadeResult,
    LayerScores,
)
from .crosstalk import (
    CrosstalkDetector,
    CrosstalkResult,
    InterferenceEdge,
)
from .drift import (
    DriftDetector,
    DriftMeasurement,
    DriftResult,
)
from .feedback import (
    FeedbackMonitor,
    FeedbackResult,
    FeedbackStatus,
)
from .gates import (
    GateCheck,
    ReleaseDecision,
    ResolutionGate,
    ResolutionResult,
)
from .stress import (
    StressFingerprinter,
    StressResult,
    StressSignature,
)
from .temporal import (
    PathwayTrajectory,
    TemporalResult,
    TemporalTracker,
    TrajectoryType,
)
from .tension import (
    TensionResult,
    TensionScorer,
)

__all__ = [
    # F1: Cascade Detection
    "CascadeAnalyzer", "CascadeResult", "LayerScores",
    # F2: Temporal Tracking
    "TemporalTracker", "TemporalResult", "PathwayTrajectory", "TrajectoryType",
    # F3: Tension Scoring
    "TensionScorer", "TensionResult",
    # F4: Resolution Gates
    "ResolutionGate", "ResolutionResult", "ReleaseDecision", "GateCheck",
    # F5: Drift Detection
    "DriftDetector", "DriftResult", "DriftMeasurement",
    # F6: Off-Target Detection
    "OffTargetDetector", "OffTargetResult", "ActivationClass", "PathwayActivationResult",
    # F7: Batch Heterogeneity
    "HeterogeneityProfiler", "HeterogeneityResult", "SubpopulationInfo",
    # F8: Dosage & Stoichiometry
    "DosageAnalyzer", "DosageAnalysisResult", "DosageState", "PathwayDosageResult", "StoichiometryResult",
    # F9: Crosstalk Detection
    "CrosstalkDetector", "CrosstalkResult", "InterferenceEdge",
    # F10: Feedback Monitoring
    "FeedbackMonitor", "FeedbackResult", "FeedbackStatus",
    # F11: Stress Fingerprinting
    "StressFingerprinter", "StressResult", "StressSignature",
    # F12: Atlas Comparison
    "AtlasComparator", "AtlasComparisonResult", "AtlasEntry", "CellMapping",
]
