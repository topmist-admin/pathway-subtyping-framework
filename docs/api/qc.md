# Molecular QC Layer API

> **Module**: `pathway_subtyping.qc`
> **Install**: `pip install pathway-subtyping[qc]`

Provides pathway-level quality control for manufactured and engineered cells, moving from "pathway activity profiling" to "pathway resolution monitoring." 12 features across two workstreams.

---

## Quick Example

```python
from pathway_subtyping.qc import (
    OffTargetDetector, DosageAnalyzer, HeterogeneityProfiler,
    CascadeAnalyzer, TensionScorer, ResolutionGate,
)

# F6: Off-target detection
detector = OffTargetDetector(
    intended_pathways=["HALLMARK_INFLAMMATORY_RESPONSE"],
    excluded_pathways=["HALLMARK_APOPTOSIS"],
)
offtarget_result = detector.scan(pathway_scores)

# F8: Dosage gating
analyzer = DosageAnalyzer(
    therapeutic_windows={"HALLMARK_INFLAMMATORY_RESPONSE": (0.4, 0.8)},
)
dosage_result = analyzer.analyze(pathway_scores)

# F4: Resolution gate (integrates all signals)
gate = ResolutionGate(min_cascade_ccr=0.7, max_heterogeneity=0.15)
decision = gate.evaluate(
    dosage_result=dosage_result,
    offtarget_result=offtarget_result,
)
print(decision.decision)  # RELEASE / HOLD / REJECT
```

---

## Feature Summary

| # | Feature | Class | KG Required | Description |
|---|---------|-------|-------------|-------------|
| F1 | Cascade Detection | `CascadeAnalyzer` | Yes | Detects incomplete signaling cascades via KG topology |
| F2 | Temporal Tracking | `TemporalTracker` | No | Trajectory classification (resolving/stalled/reversing/oscillating) |
| F3 | Tension Scoring | `TensionScorer` | No | Quantifies molecular tension from open signaling loops |
| F4 | Resolution Gates | `ResolutionGate` | No | Unified RELEASE/HOLD/REJECT decision integrating all QC signals |
| F5 | Drift Detection | `DriftDetector` | No | Cumulative pathway drift from baseline across passages |
| F6 | Off-Target Detection | `OffTargetDetector` | No | Classifies activations as INTENDED/TOLERATED/OFF_TARGET/EXCLUDED_VIOLATION |
| F7 | Heterogeneity | `HeterogeneityProfiler` | No | Batch uniformity: conformity scores, bimodality, subpopulation detection |
| F8 | Dosage & Stoichiometry | `DosageAnalyzer` | No | UNDER/IN_RANGE/OVER classification with therapeutic windows |
| F9 | Crosstalk Detection | `CrosstalkDetector` | Yes | ⚠️ **EXPERIMENTAL — do not interpret results.** See [F9 notice](#f9-experimental--do-not-interpret) |
| F10 | Feedback Monitoring | `FeedbackMonitor` | No | Activator-inhibitor correlation (intact/decoupled/inverted) |
| F11 | Stress Fingerprinting | `StressFingerprinter` | No | Matches pathway patterns to 6 known stressor signatures |
| F12 | Atlas Comparison | `AtlasComparator` | No | Distance from reference atlas with nearest-type mapping |

---

## Workstream A — KG-Independent

### `OffTargetDetector` (F6)

Detects unintended pathway activations by comparing against an engineering specification.

```python
detector = OffTargetDetector(
    intended_pathways: List[str],       # Pathways that SHOULD activate
    excluded_pathways: List[str],       # Pathways that must NOT activate (safety)
    noise_threshold: float = 0.5,       # Z-score above which activation is real
    cell_fraction_threshold: float = 0.1,  # Min fraction of cells to flag
)

result = detector.scan(pathway_scores: pd.DataFrame) -> OffTargetResult
result.n_off_target           # Count of off-target pathways
result.n_excluded_violations  # Count of safety violations
result.safety_flag            # True if excluded pathway activated
result.get_off_target_pathways()      # List of off-target pathway names
result.get_excluded_violations()      # List of violated excluded pathways
```

**Classification:** `ActivationClass.INTENDED` | `TOLERATED` | `OFF_TARGET` | `EXCLUDED_VIOLATION`

### `HeterogeneityProfiler` (F7)

Profiles batch uniformity against a target pathway profile.

```python
profiler = HeterogeneityProfiler(
    target_profile: Optional[Dict[str, float]],  # Target per-pathway scores
    conformity_threshold: float = 2.0,            # Max distance to be "conforming"
    heterogeneity_threshold: float = 0.15,        # Max non-conforming fraction
)

result = profiler.profile(pathway_scores) -> HeterogeneityResult
result.heterogeneity_index    # 0.0 (uniform) to 1.0 (heterogeneous)
result.bimodality_coefficient # >0.555 suggests bimodal distribution
result.conformity_scores      # Per-cell distance from target (numpy array)
result.subpopulations         # List[SubpopulationInfo] from DBSCAN
```

### `DosageAnalyzer` (F8)

Classifies pathway activation into three states relative to therapeutic windows.

```python
analyzer = DosageAnalyzer(
    therapeutic_windows: Dict[str, Tuple[float, float]],  # pathway -> (min, max)
    pathway_ratios: List[Dict],     # Stoichiometric ratio checks
    toxicity_weights: Dict[str, float],  # Per-pathway toxicity weight
)

result = analyzer.analyze(pathway_scores) -> DosageAnalysisResult
result.n_under / result.n_in_range / result.n_over
result.toxicity_risk_score    # Weighted over-activation risk
result.stoichiometry_results  # List[StoichiometryResult]
```

**States:** `DosageState.UNDER` | `IN_RANGE` | `OVER`

### `AtlasComparator` (F12)

Compares batch against a reference atlas (published atlas, gold-standard batch, or healthy tissue reference).

```python
atlas = AtlasComparator(distance_threshold=3.0)
atlas.add_reference("healthy_T_cell", {"PW_A": 0.6, "PW_B": 0.3})
atlas.load_atlas_csv("atlas.csv", name_column="cell_type")

result = atlas.compare(pathway_scores, expected_type="healthy_T_cell")
result.median_distance / result.p90_distance
result.fraction_on_target     # Cells matching expected type
result.type_distribution      # Dict[type_name, count]
result.cell_mappings          # Per-cell nearest type + distance
```

---

## Workstream B — KG-Dependent

### `CascadeAnalyzer` (F1)

Detects incomplete signaling cascades using KG topology to partition genes into upstream/intermediate/downstream layers.

```python
from pathway_subtyping.knowledge_graph import KnowledgeGraph

analyzer = CascadeAnalyzer(
    kg: KnowledgeGraph,
    activation_threshold: float = 0.3,
    completion_threshold: float = 0.5,
)

result = analyzer.analyze(expression, pathways=["MAPK"]) -> CascadeResult
result.mean_completion_ratio  # 0.0 (all stalled) to 1.0 (all complete)
result.n_incomplete           # Pathways with upstream > threshold but downstream < threshold
result.per_cell_ccr           # DataFrame (cells x pathways) of completion ratios
result.get_incomplete_pathways()  # List of stalled pathway names
```

### `TemporalTracker` (F2)

Tracks pathway trajectories across manufacturing timepoints.

```python
tracker = TemporalTracker(
    expected_profile: Dict[str, float],  # Target "done" state
    stall_threshold: float = 0.05,
)
tracker.register_timepoint(scores_day0, "day_0")
tracker.register_timepoint(scores_day3, "day_3")
tracker.register_timepoint(scores_day7, "day_7")

result = tracker.compute_all_trajectories() -> TemporalResult
result.n_stalled / result.n_reversing / result.n_oscillating
result.get_problem_pathways()
```

**Trajectory types:** `TrajectoryType.RESOLVING` | `STALLED` | `REVERSING` | `OSCILLATING`

### `TensionScorer` (F3)

Quantifies molecular tension: `T(cell) = sum(w * (1 - CCR) * A)`.

```python
scorer = TensionScorer(
    pathway_weights: Dict[str, float],  # Higher = more tension when stalled
    tension_threshold: float = 1.0,
)
result = scorer.score_batch(cascade_result) -> TensionResult
result.mean_tension / result.max_tension
result.per_pathway_contribution  # Which pathways contribute most
result.gate_passed               # True if mean tension <= threshold
```

### `ResolutionGate` (F4)

Unified manufacturing release gate integrating all QC signals.

```python
gate = ResolutionGate(
    min_cascade_ccr=0.7, max_tension=1.0,
    max_heterogeneity=0.15, max_offtarget=0,
)
result = gate.evaluate(
    cascade_result=..., tension_result=...,
    heterogeneity_result=..., dosage_result=..., offtarget_result=...,
)
result.decision  # ReleaseDecision.RELEASE / HOLD / REJECT
result.gate_checks  # List[GateCheck] with per-gate pass/fail
```

### `DriftDetector` (F5)

Detects cumulative pathway drift from baseline.

```python
detector = DriftDetector(max_drift=0.5)
detector.establish_baseline({"PW_A": 0.9, "PW_B": 0.8})
detector.measure_drift({"PW_A": 0.7, "PW_B": 0.6}, "passage_5")

result = detector.drift_trajectory() -> DriftResult
result.alarm_triggered / result.alarm_timepoint
result.drift_drivers  # Top pathways driving drift
```

### `CrosstalkDetector` (F9), `FeedbackMonitor` (F10), `StressFingerprinter` (F11)

#### F9: EXPERIMENTAL — do not interpret

> ⚠️ **`competition_score` does not currently measure competition at shared nodes.**
> `_compute_competition` subtracts the shared gene set from both pathways and then
> returns `mean(A-exclusive) - mean(B-exclusive)`; the shared genes never enter the
> arithmetic and act only as a presence gate.
>
> Demonstrated: varying shared-node expression across four orders of magnitude
> (0.1 → 500) leaves the score bit-identical, while holding the shared nodes
> perfectly balanced and varying only the exclusive genes swings it from +8.0 to
> −8.0 and flips `dominant`.
>
> Consequently **`competition_score`, `dominant`, and the PASS/FAIL summary are all
> uninterpretable**, and `n_significant` (which thresholds `abs(score)` at 0.3)
> would fail most real batches for an unrelated reason.
>
> F9 is therefore **soft-deprecated**: constructing the detector emits a
> `FutureWarning`, and the names are withheld from `pathway_subtyping.qc.__all__`.
> They remain importable, so existing code does not break. The replacement is
> specified in [`../roadmap-f9-competition-model.md`](../roadmap-f9-competition-model.md)
> — a partial-correlation screen confirmed by a starvation *interaction* model,
> targeted at **v0.9.0**.
>
> This does **not** affect `KnowledgeGraph.get_pathway_crosstalk()` or
> `get_shared_genes()`, which are unrelated topology helpers — nor **F9 of the v0.6
> release** (`qc.offtarget_sequence`, Evo 2 off-target scoring), a different feature
> that happens to share the label.

```python
# F9: Crosstalk — EXPERIMENTAL, emits FutureWarning; results are not interpretable
detector = CrosstalkDetector(kg, competition_threshold=0.3)
result = detector.detect(expression, pathways=["PW_A", "PW_B"])

# F10: Feedback integrity
monitor = FeedbackMonitor(loops=[("ACTIVATOR_PW", "INHIBITOR_PW")])
result = monitor.check(pathway_scores)
# FeedbackStatus.INTACT / DECOUPLED / INVERTED

# F11: Stress fingerprinting
fp = StressFingerprinter()
fp.load_default_signatures()  # 6 stressors: hypoxia, nutrient, oxidative, mechanical, thermal, pH
result = fp.fingerprint(pathway_scores)
result.primary_stressor  # e.g., "hypoxia"
result.remediation        # e.g., "Check dissolved oxygen in bioreactor"
```

---

## Test Infrastructure

The QC layer includes a full synthetic manufacturing simulator for validation:

```python
from pathway_subtyping.qc.testing import ManufacturingSimulator, ManufacturingSpec

sim = ManufacturingSimulator(seed=42)
batch = sim.generate_healthy_batch(n_cells=500, spec=spec)
injected = sim.inject_offtarget(batch, pathways=["HALLMARK_APOPTOSIS"], cell_fraction=0.2)
```

9 injectable defect types, 12x12 orthogonality matrix, severity titration harness, and 3 end-to-end scenario tests (CAR T, neural organoid, 20-passage stability).

### Multi-Level Validation (L2/L3/L4)

Beyond synthetic testing, the QC layer includes a 3-level validation framework:

| Level | Module | Purpose |
|-------|--------|---------|
| L2 Retrospective | `qc.testing.retrospective` | Correlate QC scores with clinical outcomes on public datasets (TCGA-COAD, GSE65682, GSE15402) |
| L3 Scenarios | `qc.testing.scenarios` | End-to-end scenario execution with defect injection, detection validation, and ground truth checking |
| L4 Prospective | `qc.testing.prospective` | Shadow protocol for non-interventional PSF alongside standard QC, with concordance analysis and threshold calibration |

---

Research use only. Not for clinical decision-making.
