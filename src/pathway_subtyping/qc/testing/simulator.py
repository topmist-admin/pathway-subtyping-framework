"""
Manufacturing Simulator for Molecular QC Testing.

Generates synthetic expression data with controllable, labeled defects
so every QC feature can be tested with exact sensitivity/specificity.

Defect types:
    - stalled_cascade: Upstream active, downstream suppressed
    - temporal_drift: Progressive CCR degradation across passages
    - offtarget: Unintended pathway activation
    - overdose: Pathway activation above therapeutic window
    - subpopulation: Non-conforming cell subpopulation
    - crosstalk: Shared node competition between pathways
    - broken_feedback: Decoupled or inverted activator-inhibitor loops
    - stress: Environmental stress signature overlay
    - atlas_drift: Shift from reference atlas centroid

Research use only. Not for clinical decision-making.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class DefectType(Enum):
    """Types of injectable manufacturing defects."""

    STALLED_CASCADE = "stalled_cascade"
    TEMPORAL_DRIFT = "temporal_drift"
    OFFTARGET = "offtarget"
    OVERDOSE = "overdose"
    SUBPOPULATION = "subpopulation"
    CROSSTALK = "crosstalk"
    BROKEN_FEEDBACK = "broken_feedback"
    STRESS = "stress"
    ATLAS_DRIFT = "atlas_drift"


class StressorType(Enum):
    """Types of environmental stressors."""

    HYPOXIA = "hypoxia"
    NUTRIENT_DEPLETION = "nutrient_depletion"
    OXIDATIVE = "oxidative"
    MECHANICAL = "mechanical"
    THERMAL = "thermal"
    PH = "ph"


class FeedbackBreakType(Enum):
    """Types of feedback loop failures."""

    DECOUPLED = "decoupled"
    INVERTED = "inverted"


class CompetitionDirection(Enum):
    """Direction of pathway crosstalk competition."""

    A_DOMINATES = "a_dominates"
    B_DOMINATES = "b_dominates"
    BALANCED = "balanced"


@dataclass
class DefectLabel:
    """Per-cell defect label with metadata."""

    affected: bool
    defect_type: Optional[DefectType] = None
    severity: float = 0.0
    affected_pathways: List[str] = field(default_factory=list)


@dataclass
class DefectMetadata:
    """Metadata describing an injected defect."""

    defect_type: DefectType
    severity: float
    affected_pathways: List[str]
    parameters: Dict[str, Any] = field(default_factory=dict)
    expected_detections: Dict[str, bool] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "defect_type": self.defect_type.value,
            "severity": self.severity,
            "affected_pathways": self.affected_pathways,
            "parameters": self.parameters,
            "expected_detections": self.expected_detections,
        }


@dataclass
class ManufacturingSpec:
    """Engineering specification for a manufactured cell product.

    Attributes:
        intended_pathways: Pathways that SHOULD activate
        excluded_pathways: Pathways that must NOT activate (safety-critical)
        therapeutic_windows: Per-pathway (min, max) activation ranges
        pathway_ratios: Expected stoichiometric ratios between pathway pairs
        feedback_loops: Activator-inhibitor pairs that must remain coupled
        target_profile: Target pathway activation profile for conformity scoring
    """

    intended_pathways: List[str] = field(default_factory=list)
    excluded_pathways: List[str] = field(default_factory=list)
    therapeutic_windows: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    pathway_ratios: Dict[Tuple[str, str], Tuple[float, float]] = field(
        default_factory=dict
    )
    feedback_loops: List[Tuple[str, str]] = field(default_factory=list)
    target_profile: Optional[Dict[str, float]] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "intended_pathways": self.intended_pathways,
            "excluded_pathways": self.excluded_pathways,
            "therapeutic_windows": {
                k: list(v) for k, v in self.therapeutic_windows.items()
            },
            "feedback_loops": [list(pair) for pair in self.feedback_loops],
        }


@dataclass
class InjectedBatch:
    """A batch of synthetic expression data with injected defects.

    Attributes:
        expression: Expression matrix (cells x genes)
        pathway_scores: Pathway activation scores (cells x pathways)
        cell_labels: Per-cell defect labels
        defect_metadata: Description of what was injected
        spec: The manufacturing specification used
    """

    expression: pd.DataFrame
    pathway_scores: pd.DataFrame
    cell_labels: List[DefectLabel]
    defect_metadata: List[DefectMetadata]
    spec: ManufacturingSpec

    @property
    def n_cells(self) -> int:
        return len(self.expression)

    @property
    def n_genes(self) -> int:
        return self.expression.shape[1]

    @property
    def n_pathways(self) -> int:
        return self.pathway_scores.shape[1]

    @property
    def affected_mask(self) -> np.ndarray:
        """Boolean array: True for cells affected by any defect."""
        return np.array([label.affected for label in self.cell_labels])

    def to_dict(self) -> Dict[str, Any]:
        return {
            "n_cells": self.n_cells,
            "n_genes": self.n_genes,
            "n_pathways": self.n_pathways,
            "n_affected": int(self.affected_mask.sum()),
            "defects": [m.to_dict() for m in self.defect_metadata],
        }


class ManufacturingSimulator:
    """Generates synthetic manufacturing batches with injectable defects.

    Usage::

        sim = ManufacturingSimulator(seed=42)
        spec = ManufacturingSpec(
            intended_pathways=["HALLMARK_INFLAMMATORY_RESPONSE"],
            excluded_pathways=["HALLMARK_APOPTOSIS"],
            therapeutic_windows={"HALLMARK_INFLAMMATORY_RESPONSE": (0.4, 0.8)},
        )
        batch = sim.generate_healthy_batch(n_cells=500, n_genes=200, spec=spec)
        injected = sim.inject_offtarget(batch, pathways=["HALLMARK_APOPTOSIS"],
                                         cell_fraction=0.2, activation_level=0.7)
    """

    # Default Hallmark pathway names for synthetic data
    DEFAULT_PATHWAYS = [
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_COMPLEMENT",
        "HALLMARK_APOPTOSIS",
        "HALLMARK_P53_PATHWAY",
        "HALLMARK_HYPOXIA",
        "HALLMARK_GLYCOLYSIS",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "HALLMARK_FATTY_ACID_METABOLISM",
        "HALLMARK_MTORC1_SIGNALING",
        "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
        "HALLMARK_KRAS_SIGNALING_UP",
        "HALLMARK_KRAS_SIGNALING_DN",
        "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
        "HALLMARK_NOTCH_SIGNALING",
        "HALLMARK_HEDGEHOG_SIGNALING",
        "HALLMARK_TGF_BETA_SIGNALING",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_ANGIOGENESIS",
        "HALLMARK_MYOGENESIS",
        "HALLMARK_ADIPOGENESIS",
        "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    ]

    # Stress signatures: activated and suppressed pathways per stressor
    STRESS_SIGNATURES: Dict[StressorType, Dict[str, List[str]]] = {
        StressorType.HYPOXIA: {
            "activated": ["HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "HALLMARK_ANGIOGENESIS"],
            "suppressed": ["HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM"],
        },
        StressorType.NUTRIENT_DEPLETION: {
            "activated": ["HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_P53_PATHWAY"],
            "suppressed": ["HALLMARK_MTORC1_SIGNALING", "HALLMARK_ADIPOGENESIS"],
        },
        StressorType.OXIDATIVE: {
            "activated": ["HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS"],
            "suppressed": ["HALLMARK_OXIDATIVE_PHOSPHORYLATION"],
        },
        StressorType.MECHANICAL: {
            "activated": ["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_MYOGENESIS"],
            "suppressed": ["HALLMARK_ADIPOGENESIS"],
        },
        StressorType.THERMAL: {
            "activated": ["HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_P53_PATHWAY"],
            "suppressed": ["HALLMARK_MTORC1_SIGNALING"],
        },
        StressorType.PH: {
            "activated": ["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"],
            "suppressed": ["HALLMARK_MTORC1_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"],
        },
    }

    def __init__(
        self,
        pathways: Optional[List[str]] = None,
        genes_per_pathway: int = 20,
        seed: Optional[int] = None,
    ):
        self.pathways = pathways or self.DEFAULT_PATHWAYS
        self.genes_per_pathway = genes_per_pathway
        self.seed = seed
        self._rng = np.random.default_rng(seed)

        # Build gene-pathway mapping
        self._pathway_genes: Dict[str, List[str]] = {}
        gene_idx = 0
        for pw in self.pathways:
            genes = [f"GENE_{gene_idx + i}" for i in range(genes_per_pathway)]
            self._pathway_genes[pw] = genes
            gene_idx += genes_per_pathway
        self._all_genes = [
            g for genes in self._pathway_genes.values() for g in genes
        ]

    @property
    def pathway_genes(self) -> Dict[str, List[str]]:
        """Gene-pathway mapping used by this simulator."""
        return dict(self._pathway_genes)

    @property
    def all_genes(self) -> List[str]:
        """All gene names in the simulation."""
        return list(self._all_genes)

    def generate_healthy_batch(
        self,
        n_cells: int = 500,
        n_genes: Optional[int] = None,
        spec: Optional[ManufacturingSpec] = None,
        base_expression: float = 5.0,
        noise_std: float = 1.0,
    ) -> InjectedBatch:
        """Generate a healthy batch where all cascades are complete and on-target.

        Args:
            n_cells: Number of cells in the batch.
            n_genes: Ignored (gene count determined by pathway structure).
            spec: Manufacturing specification. If None, a default spec is created.
            base_expression: Mean expression level for all genes.
            noise_std: Standard deviation of expression noise.

        Returns:
            InjectedBatch with no defects injected.
        """
        if spec is None:
            spec = ManufacturingSpec(
                intended_pathways=self.pathways[:5],
                excluded_pathways=self.pathways[5:8],
            )

        logger.info(
            "[QC Simulator] Generating healthy batch: %d cells, %d genes, %d pathways",
            n_cells,
            len(self._all_genes),
            len(self.pathways),
        )

        # Generate expression: base + noise, with intended pathways slightly elevated
        expression = np.full(
            (n_cells, len(self._all_genes)), base_expression, dtype=np.float64
        )
        expression += self._rng.normal(0, noise_std, size=expression.shape)

        # Elevate intended pathways
        for pw in spec.intended_pathways:
            if pw in self._pathway_genes:
                col_indices = [
                    self._all_genes.index(g) for g in self._pathway_genes[pw]
                ]
                expression[:, col_indices] += 2.0

        # Suppress excluded pathways
        for pw in spec.excluded_pathways:
            if pw in self._pathway_genes:
                col_indices = [
                    self._all_genes.index(g) for g in self._pathway_genes[pw]
                ]
                expression[:, col_indices] -= 2.0

        expr_df = pd.DataFrame(
            expression,
            index=[f"cell_{i}" for i in range(n_cells)],
            columns=self._all_genes,
        )

        # Compute pathway scores as mean expression per pathway
        pathway_scores = self._compute_pathway_scores(expr_df)

        # All cells healthy
        cell_labels = [DefectLabel(affected=False) for _ in range(n_cells)]

        return InjectedBatch(
            expression=expr_df,
            pathway_scores=pathway_scores,
            cell_labels=cell_labels,
            defect_metadata=[],
            spec=spec,
        )

    def inject_stalled_cascades(
        self,
        batch: InjectedBatch,
        pathways: List[str],
        layer: str = "downstream",
        severity: float = 0.5,
    ) -> InjectedBatch:
        """Inject stalled cascades: upstream active, downstream suppressed.

        Args:
            batch: Source batch to inject into.
            pathways: Pathways to stall.
            layer: Which layer to suppress ("downstream" or "intermediate").
            severity: 0.0 (subtle) to 1.0 (complete stall).

        Returns:
            New InjectedBatch with stalled cascade defect.
        """
        expr = batch.expression.values.copy()
        n_cells = len(expr)
        affected = np.ones(n_cells, dtype=bool)

        for pw in pathways:
            if pw not in self._pathway_genes:
                continue
            genes = self._pathway_genes[pw]
            n_genes = len(genes)
            col_indices = [self._all_genes.index(g) for g in genes]

            # Split into upstream (first third), middle, downstream (last third)
            n_up = max(1, n_genes // 3)
            n_down = max(1, n_genes // 3)

            if layer == "downstream":
                suppress_indices = col_indices[-n_down:]
            else:
                suppress_indices = col_indices[n_up : n_up + (n_genes - n_up - n_down)]

            # Elevate upstream, suppress target layer
            upstream_indices = col_indices[:n_up]
            expr[:, upstream_indices] += 2.0 * severity
            expr[:, suppress_indices] -= 3.0 * severity

        expr_df = pd.DataFrame(expr, index=batch.expression.index, columns=batch.expression.columns)
        pathway_scores = self._compute_pathway_scores(expr_df)

        cell_labels = [
            DefectLabel(
                affected=True,
                defect_type=DefectType.STALLED_CASCADE,
                severity=severity,
                affected_pathways=pathways,
            )
            for _ in range(n_cells)
        ]

        metadata = DefectMetadata(
            defect_type=DefectType.STALLED_CASCADE,
            severity=severity,
            affected_pathways=pathways,
            parameters={"layer": layer},
            expected_detections={
                "F1_cascade": True,
                "F3_tension": True,
                "F4_resolution_gate": True,
            },
        )

        return InjectedBatch(
            expression=expr_df,
            pathway_scores=pathway_scores,
            cell_labels=cell_labels,
            defect_metadata=batch.defect_metadata + [metadata],
            spec=batch.spec,
        )

    def inject_temporal_drift(
        self,
        batch: InjectedBatch,
        n_passages: int = 10,
        drift_rate: float = 0.02,
        affected_pathways: Optional[List[str]] = None,
    ) -> List[InjectedBatch]:
        """Inject progressive pathway drift across passages.

        Args:
            batch: Baseline batch (passage 0).
            n_passages: Number of passages to simulate.
            drift_rate: CCR degradation per passage per affected pathway.
            affected_pathways: Pathways that drift. Defaults to first 2 intended.

        Returns:
            List of InjectedBatch (one per passage, including baseline).
        """
        if affected_pathways is None:
            affected_pathways = batch.spec.intended_pathways[:2]

        passages = [batch]
        base_expr = batch.expression.values.copy()

        for p in range(1, n_passages + 1):
            expr = base_expr.copy()
            drift_amount = drift_rate * p

            for pw in affected_pathways:
                if pw not in self._pathway_genes:
                    continue
                col_indices = [
                    self._all_genes.index(g) for g in self._pathway_genes[pw]
                ]
                expr[:, col_indices] -= drift_amount * 5.0
                expr[:, col_indices] += self._rng.normal(
                    0, 0.1 * p, size=(len(expr), len(col_indices))
                )

            expr_df = pd.DataFrame(
                expr, index=batch.expression.index, columns=batch.expression.columns
            )
            pathway_scores = self._compute_pathway_scores(expr_df)

            severity = min(1.0, drift_amount / 0.2)
            cell_labels = [
                DefectLabel(
                    affected=p > 0,
                    defect_type=DefectType.TEMPORAL_DRIFT,
                    severity=severity,
                    affected_pathways=affected_pathways,
                )
                for _ in range(len(expr))
            ]

            metadata = DefectMetadata(
                defect_type=DefectType.TEMPORAL_DRIFT,
                severity=severity,
                affected_pathways=affected_pathways,
                parameters={"passage": p, "drift_rate": drift_rate},
                expected_detections={
                    "F2_temporal": True,
                    "F5_drift": True,
                    "F12_atlas": severity > 0.3,
                },
            )

            passages.append(
                InjectedBatch(
                    expression=expr_df,
                    pathway_scores=pathway_scores,
                    cell_labels=cell_labels,
                    defect_metadata=[metadata],
                    spec=batch.spec,
                )
            )

        logger.info(
            "[QC Simulator] Injected temporal drift: %d passages, rate=%.3f, %d pathways",
            n_passages,
            drift_rate,
            len(affected_pathways),
        )
        return passages

    def inject_offtarget(
        self,
        batch: InjectedBatch,
        pathways: List[str],
        cell_fraction: float = 0.2,
        activation_level: float = 0.7,
    ) -> InjectedBatch:
        """Inject off-target pathway activation in a fraction of cells.

        Args:
            batch: Source batch.
            pathways: Off-target pathways to activate (should NOT be in spec).
            cell_fraction: Fraction of cells showing off-target activation.
            activation_level: Strength of activation (0-1 scale).

        Returns:
            New InjectedBatch with off-target defect.
        """
        expr = batch.expression.values.copy()
        n_cells = len(expr)
        n_affected = max(1, int(n_cells * cell_fraction))

        affected_indices = self._rng.choice(n_cells, size=n_affected, replace=False)
        affected_mask = np.zeros(n_cells, dtype=bool)
        affected_mask[affected_indices] = True

        for pw in pathways:
            if pw not in self._pathway_genes:
                continue
            col_indices = [
                self._all_genes.index(g) for g in self._pathway_genes[pw]
            ]
            expr[affected_indices[:, None], col_indices] += 4.0 * activation_level

        expr_df = pd.DataFrame(
            expr, index=batch.expression.index, columns=batch.expression.columns
        )
        pathway_scores = self._compute_pathway_scores(expr_df)

        cell_labels = [
            DefectLabel(
                affected=bool(affected_mask[i]),
                defect_type=DefectType.OFFTARGET if affected_mask[i] else None,
                severity=activation_level if affected_mask[i] else 0.0,
                affected_pathways=pathways if affected_mask[i] else [],
            )
            for i in range(n_cells)
        ]

        metadata = DefectMetadata(
            defect_type=DefectType.OFFTARGET,
            severity=activation_level,
            affected_pathways=pathways,
            parameters={"cell_fraction": cell_fraction},
            expected_detections={
                "F6_offtarget": True,
            },
        )

        logger.info(
            "[QC Simulator] Injected off-target: %d/%d cells, %d pathways, level=%.2f",
            n_affected,
            n_cells,
            len(pathways),
            activation_level,
        )

        return InjectedBatch(
            expression=expr_df,
            pathway_scores=pathway_scores,
            cell_labels=cell_labels,
            defect_metadata=batch.defect_metadata + [metadata],
            spec=batch.spec,
        )

    def inject_overdose(
        self,
        batch: InjectedBatch,
        pathways: List[str],
        multiplier: float = 2.0,
    ) -> InjectedBatch:
        """Inject over-activation of pathways above therapeutic window.

        Args:
            batch: Source batch.
            pathways: Pathways to over-activate.
            multiplier: Scale factor above normal (1.5x, 2x, 5x).

        Returns:
            New InjectedBatch with overdose defect.
        """
        expr = batch.expression.values.copy()
        n_cells = len(expr)

        for pw in pathways:
            if pw not in self._pathway_genes:
                continue
            col_indices = [
                self._all_genes.index(g) for g in self._pathway_genes[pw]
            ]
            # Scale up expression for these pathways
            current = expr[:, col_indices]
            expr[:, col_indices] = current * multiplier

        expr_df = pd.DataFrame(
            expr, index=batch.expression.index, columns=batch.expression.columns
        )
        pathway_scores = self._compute_pathway_scores(expr_df)

        severity = min(1.0, (multiplier - 1.0) / 4.0)
        cell_labels = [
            DefectLabel(
                affected=True,
                defect_type=DefectType.OVERDOSE,
                severity=severity,
                affected_pathways=pathways,
            )
            for _ in range(n_cells)
        ]

        metadata = DefectMetadata(
            defect_type=DefectType.OVERDOSE,
            severity=severity,
            affected_pathways=pathways,
            parameters={"multiplier": multiplier},
            expected_detections={
                "F8_dosage": True,
                "F3_tension": multiplier > 2.0,
            },
        )

        logger.info(
            "[QC Simulator] Injected overdose: %.1fx on %d pathways",
            multiplier,
            len(pathways),
        )

        return InjectedBatch(
            expression=expr_df,
            pathway_scores=pathway_scores,
            cell_labels=cell_labels,
            defect_metadata=batch.defect_metadata + [metadata],
            spec=batch.spec,
        )

    def inject_subpopulation(
        self,
        batch: InjectedBatch,
        fraction: float = 0.1,
        deviant_profile: Optional[Dict[str, float]] = None,
    ) -> InjectedBatch:
        """Inject a non-conforming subpopulation.

        Args:
            batch: Source batch.
            fraction: Fraction of cells to replace with deviant profile.
            deviant_profile: Per-pathway deviation amounts. Defaults to +-2.0
                on first 3 pathways.

        Returns:
            New InjectedBatch with subpopulation defect.
        """
        if deviant_profile is None:
            deviant_profile = {}
            for i, pw in enumerate(self.pathways[:3]):
                deviant_profile[pw] = 3.0 if i % 2 == 0 else -3.0

        expr = batch.expression.values.copy()
        n_cells = len(expr)
        n_deviant = max(1, int(n_cells * fraction))

        deviant_indices = self._rng.choice(n_cells, size=n_deviant, replace=False)
        deviant_mask = np.zeros(n_cells, dtype=bool)
        deviant_mask[deviant_indices] = True

        for pw, shift in deviant_profile.items():
            if pw not in self._pathway_genes:
                continue
            col_indices = [
                self._all_genes.index(g) for g in self._pathway_genes[pw]
            ]
            expr[deviant_indices[:, None], col_indices] += shift

        expr_df = pd.DataFrame(
            expr, index=batch.expression.index, columns=batch.expression.columns
        )
        pathway_scores = self._compute_pathway_scores(expr_df)

        cell_labels = [
            DefectLabel(
                affected=bool(deviant_mask[i]),
                defect_type=DefectType.SUBPOPULATION if deviant_mask[i] else None,
                severity=fraction if deviant_mask[i] else 0.0,
                affected_pathways=list(deviant_profile.keys()) if deviant_mask[i] else [],
            )
            for i in range(n_cells)
        ]

        metadata = DefectMetadata(
            defect_type=DefectType.SUBPOPULATION,
            severity=fraction,
            affected_pathways=list(deviant_profile.keys()),
            parameters={"fraction": fraction, "deviant_profile": deviant_profile},
            expected_detections={
                "F7_heterogeneity": True,
                "F12_atlas": fraction > 0.15,
            },
        )

        logger.info(
            "[QC Simulator] Injected subpopulation: %d/%d deviant cells",
            n_deviant,
            n_cells,
        )

        return InjectedBatch(
            expression=expr_df,
            pathway_scores=pathway_scores,
            cell_labels=cell_labels,
            defect_metadata=batch.defect_metadata + [metadata],
            spec=batch.spec,
        )

    def inject_crosstalk(
        self,
        batch: InjectedBatch,
        pathway_a: str,
        pathway_b: str,
        shared_nodes: Optional[List[str]] = None,
        competition_direction: CompetitionDirection = CompetitionDirection.A_DOMINATES,
    ) -> InjectedBatch:
        """Inject pathway crosstalk via shared node competition.

        Creates synthetic shared genes between two pathways and skews their
        expression toward one pathway's downstream targets.

        Args:
            batch: Source batch.
            pathway_a: First competing pathway.
            pathway_b: Second competing pathway.
            shared_nodes: Explicit shared gene names. If None, creates synthetic ones.
            competition_direction: Which pathway dominates the shared nodes.

        Returns:
            New InjectedBatch with crosstalk defect.
        """
        expr = batch.expression.values.copy()
        n_cells = len(expr)

        genes_a = self._pathway_genes.get(pathway_a, [])
        genes_b = self._pathway_genes.get(pathway_b, [])
        if not genes_a or not genes_b:
            logger.warning("[QC Simulator] Pathways not found for crosstalk injection")
            return batch

        # Use middle genes of each pathway as "shared nodes"
        n_shared = min(3, len(genes_a) // 3, len(genes_b) // 3)
        mid_a = len(genes_a) // 2
        mid_b = len(genes_b) // 2
        shared_a_indices = [self._all_genes.index(genes_a[mid_a + i]) for i in range(n_shared)]
        shared_b_indices = [self._all_genes.index(genes_b[mid_b + i]) for i in range(n_shared)]

        # Downstream genes (last third)
        n_down = max(1, len(genes_a) // 3)
        down_a = [self._all_genes.index(genes_a[-i - 1]) for i in range(n_down)]
        down_b = [self._all_genes.index(genes_b[-i - 1]) for i in range(n_down)]

        if competition_direction == CompetitionDirection.A_DOMINATES:
            # A's shared nodes are high, B's downstream is starved
            expr[:, shared_a_indices] += 2.0
            expr[:, down_b] -= 2.5
        elif competition_direction == CompetitionDirection.B_DOMINATES:
            expr[:, shared_b_indices] += 2.0
            expr[:, down_a] -= 2.5
        else:
            # Balanced — mild competition both ways
            expr[:, shared_a_indices] += 1.0
            expr[:, shared_b_indices] += 1.0
            expr[:, down_a] -= 0.5
            expr[:, down_b] -= 0.5

        expr_df = pd.DataFrame(
            expr, index=batch.expression.index, columns=batch.expression.columns
        )
        pathway_scores = self._compute_pathway_scores(expr_df)

        cell_labels = [
            DefectLabel(
                affected=True,
                defect_type=DefectType.CROSSTALK,
                severity=0.7,
                affected_pathways=[pathway_a, pathway_b],
            )
            for _ in range(n_cells)
        ]

        metadata = DefectMetadata(
            defect_type=DefectType.CROSSTALK,
            severity=0.7,
            affected_pathways=[pathway_a, pathway_b],
            parameters={
                "competition_direction": competition_direction.value,
                "n_shared_nodes": n_shared,
            },
            expected_detections={
                "F9_crosstalk": True,
                "F1_cascade": True,
            },
        )

        logger.info(
            "[QC Simulator] Injected crosstalk: %s vs %s (%s)",
            pathway_a,
            pathway_b,
            competition_direction.value,
        )

        return InjectedBatch(
            expression=expr_df,
            pathway_scores=pathway_scores,
            cell_labels=cell_labels,
            defect_metadata=batch.defect_metadata + [metadata],
            spec=batch.spec,
        )

    def inject_broken_feedback(
        self,
        batch: InjectedBatch,
        loops: List[Tuple[str, str]],
        break_type: FeedbackBreakType = FeedbackBreakType.DECOUPLED,
    ) -> InjectedBatch:
        """Inject broken feedback loops.

        For each (activator_pathway, inhibitor_pathway) pair, decouple or
        invert the correlation between their expression.

        Args:
            batch: Source batch.
            loops: List of (activator_pathway, inhibitor_pathway) pairs.
            break_type: DECOUPLED (inhibitor independent) or INVERTED
                (inhibitor suppressed when activator is high).

        Returns:
            New InjectedBatch with broken feedback defect.
        """
        expr = batch.expression.values.copy()
        n_cells = len(expr)
        affected_pathways = []

        for activator_pw, inhibitor_pw in loops:
            if activator_pw not in self._pathway_genes or inhibitor_pw not in self._pathway_genes:
                continue
            affected_pathways.extend([activator_pw, inhibitor_pw])

            act_cols = [
                self._all_genes.index(g) for g in self._pathway_genes[activator_pw]
            ]
            inh_cols = [
                self._all_genes.index(g) for g in self._pathway_genes[inhibitor_pw]
            ]

            # Make activator high
            expr[:, act_cols] += 2.0

            if break_type == FeedbackBreakType.DECOUPLED:
                # Inhibitor is random noise — uncorrelated with activator
                expr[:, inh_cols] = self._rng.normal(5.0, 1.5, size=(n_cells, len(inh_cols)))
            else:
                # Inverted: inhibitor goes DOWN when activator goes UP
                expr[:, inh_cols] -= 2.5

        expr_df = pd.DataFrame(
            expr, index=batch.expression.index, columns=batch.expression.columns
        )
        pathway_scores = self._compute_pathway_scores(expr_df)

        cell_labels = [
            DefectLabel(
                affected=True,
                defect_type=DefectType.BROKEN_FEEDBACK,
                severity=0.8,
                affected_pathways=affected_pathways,
            )
            for _ in range(n_cells)
        ]

        metadata = DefectMetadata(
            defect_type=DefectType.BROKEN_FEEDBACK,
            severity=0.8,
            affected_pathways=affected_pathways,
            parameters={
                "break_type": break_type.value,
                "loops": [(a, i) for a, i in loops],
            },
            expected_detections={
                "F10_feedback": True,
            },
        )

        logger.info(
            "[QC Simulator] Injected broken feedback: %d loops, type=%s",
            len(loops),
            break_type.value,
        )

        return InjectedBatch(
            expression=expr_df,
            pathway_scores=pathway_scores,
            cell_labels=cell_labels,
            defect_metadata=batch.defect_metadata + [metadata],
            spec=batch.spec,
        )

    def inject_stress(
        self,
        batch: InjectedBatch,
        stressor: StressorType,
        severity: float = 0.5,
    ) -> InjectedBatch:
        """Inject environmental stress signature.

        Overlays a stressor-specific pathway signature (activated + suppressed
        pathways) from the reference stress signature database.

        Args:
            batch: Source batch.
            stressor: Type of environmental stressor.
            severity: 0.0 (subtle) to 1.0 (extreme stress).

        Returns:
            New InjectedBatch with stress defect.
        """
        sig = self.STRESS_SIGNATURES.get(stressor)
        if sig is None:
            logger.warning("[QC Simulator] Unknown stressor: %s", stressor)
            return batch

        expr = batch.expression.values.copy()
        n_cells = len(expr)
        affected_pathways = sig["activated"] + sig["suppressed"]

        for pw in sig["activated"]:
            if pw not in self._pathway_genes:
                continue
            col_indices = [
                self._all_genes.index(g) for g in self._pathway_genes[pw]
            ]
            expr[:, col_indices] += 3.0 * severity

        for pw in sig["suppressed"]:
            if pw not in self._pathway_genes:
                continue
            col_indices = [
                self._all_genes.index(g) for g in self._pathway_genes[pw]
            ]
            expr[:, col_indices] -= 3.0 * severity

        expr_df = pd.DataFrame(
            expr, index=batch.expression.index, columns=batch.expression.columns
        )
        pathway_scores = self._compute_pathway_scores(expr_df)

        cell_labels = [
            DefectLabel(
                affected=True,
                defect_type=DefectType.STRESS,
                severity=severity,
                affected_pathways=affected_pathways,
            )
            for _ in range(n_cells)
        ]

        metadata = DefectMetadata(
            defect_type=DefectType.STRESS,
            severity=severity,
            affected_pathways=affected_pathways,
            parameters={"stressor": stressor.value},
            expected_detections={
                "F11_stress": True,
                "F3_tension": severity > 0.5,
            },
        )

        logger.info(
            "[QC Simulator] Injected stress: %s, severity=%.2f",
            stressor.value,
            severity,
        )

        return InjectedBatch(
            expression=expr_df,
            pathway_scores=pathway_scores,
            cell_labels=cell_labels,
            defect_metadata=batch.defect_metadata + [metadata],
            spec=batch.spec,
        )

    def inject_atlas_drift(
        self,
        batch: InjectedBatch,
        reference_profile: Optional[Dict[str, float]] = None,
        drift_distance: float = 0.5,
    ) -> InjectedBatch:
        """Inject drift away from a reference atlas profile.

        Args:
            batch: Source batch.
            reference_profile: Target pathway profile to drift FROM.
                If None, uses current batch mean as implicit reference.
            drift_distance: Magnitude of drift (0-1 scale).

        Returns:
            New InjectedBatch with atlas drift defect.
        """
        expr = batch.expression.values.copy()
        n_cells = len(expr)

        # Drift: shift random pathways away from their current mean
        n_drift_pathways = max(1, int(len(self.pathways) * drift_distance))
        drift_pathways = list(
            self._rng.choice(self.pathways, size=n_drift_pathways, replace=False)
        )

        for pw in drift_pathways:
            if pw not in self._pathway_genes:
                continue
            col_indices = [
                self._all_genes.index(g) for g in self._pathway_genes[pw]
            ]
            direction = self._rng.choice([-1, 1])
            expr[:, col_indices] += direction * 3.0 * drift_distance

        expr_df = pd.DataFrame(
            expr, index=batch.expression.index, columns=batch.expression.columns
        )
        pathway_scores = self._compute_pathway_scores(expr_df)

        cell_labels = [
            DefectLabel(
                affected=True,
                defect_type=DefectType.ATLAS_DRIFT,
                severity=drift_distance,
                affected_pathways=drift_pathways,
            )
            for _ in range(n_cells)
        ]

        metadata = DefectMetadata(
            defect_type=DefectType.ATLAS_DRIFT,
            severity=drift_distance,
            affected_pathways=drift_pathways,
            parameters={"drift_distance": drift_distance},
            expected_detections={
                "F12_atlas": True,
                "F7_heterogeneity": drift_distance > 0.5,
            },
        )

        logger.info(
            "[QC Simulator] Injected atlas drift: distance=%.2f, %d pathways shifted",
            drift_distance,
            n_drift_pathways,
        )

        return InjectedBatch(
            expression=expr_df,
            pathway_scores=pathway_scores,
            cell_labels=cell_labels,
            defect_metadata=batch.defect_metadata + [metadata],
            spec=batch.spec,
        )

    def inject_combined(
        self,
        batch: InjectedBatch,
        defect_list: List[Tuple[str, Dict[str, Any]]],
    ) -> InjectedBatch:
        """Apply multiple defects sequentially.

        Args:
            batch: Source batch.
            defect_list: List of (defect_method_name, kwargs) tuples.
                Example: [("inject_offtarget", {"pathways": [...], "cell_fraction": 0.2}),
                          ("inject_stress", {"stressor": StressorType.HYPOXIA})]

        Returns:
            InjectedBatch with all defects applied.
        """
        result = batch
        for method_name, kwargs in defect_list:
            method = getattr(self, method_name, None)
            if method is None:
                logger.warning("[QC Simulator] Unknown injection method: %s", method_name)
                continue
            result = method(result, **kwargs)
        return result

    def _compute_pathway_scores(self, expression: pd.DataFrame) -> pd.DataFrame:
        """Compute pathway scores as z-scored mean expression per pathway."""
        scores = {}
        for pw, genes in self._pathway_genes.items():
            present = [g for g in genes if g in expression.columns]
            if present:
                pw_expr = expression[present].values
                # Z-score per gene across cells, then mean per pathway
                means = pw_expr.mean(axis=0, keepdims=True)
                stds = pw_expr.std(axis=0, keepdims=True)
                stds[stds == 0] = 1.0
                z_scored = (pw_expr - means) / stds
                scores[pw] = z_scored.mean(axis=1)
            else:
                scores[pw] = np.zeros(len(expression))

        return pd.DataFrame(scores, index=expression.index)
