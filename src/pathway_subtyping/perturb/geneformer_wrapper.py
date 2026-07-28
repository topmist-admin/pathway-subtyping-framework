"""
Geneformer wrapper for in-silico gene perturbation.

``GeneformerPerturber`` is the high-level interface used throughout
``pathway_subtyping.perturb``. It delegates actual inference to a
pluggable ``GeneformerBackend``:

    - ``OfficialBackend`` — the production path, wraps the published
      Geneformer model (Theodoris et al. 2023). Requires the optional
      ``[perturb]`` extra (torch, transformers, geneformer package, and
      a cached checkpoint). Lazy-imports; a clear ImportError is raised
      if the extra is missing.
    - ``FallbackPerturber`` — deterministic numpy substitute that
      preserves the public API so tests and development can run without
      the optional dependencies. It models a knockout as the removal of
      that gene's contribution to a PCA embedding of the expression
      matrix, and over-expression as its amplification.

The fallback is not a substitute for the real model on biological claims.
It is a correctness and interface harness.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class PerturbationMode(str, Enum):
    """Supported perturbation modes."""

    KNOCKOUT = "knockout"
    OVEREXPRESS = "overexpress"


# --------------------------------------------------------------------------- #
# Backend ABC + official stub
# --------------------------------------------------------------------------- #


class GeneformerBackend:
    """Abstract base class for Geneformer inference backends.

    Subclasses implement ``embed`` and ``perturb``. Custom backends must
    pass the golden-test parity harness (``tests/test_perturb.py``
    ``@pytest.mark.perturb_golden``) before being used in production.
    """

    embedding_dim: int

    def embed(self, expression: pd.DataFrame) -> np.ndarray:
        raise NotImplementedError

    def perturb(
        self,
        expression: pd.DataFrame,
        gene: str,
        mode: PerturbationMode,
    ) -> np.ndarray:
        raise NotImplementedError


class OfficialBackend(GeneformerBackend):
    """Wraps the published Geneformer model (Theodoris et al. 2023).

    Production-path perturbation backend. Requires the optional
    ``[perturb]`` extra:

        pip install pathway-subtyping[perturb]

    plus a local Geneformer V2 checkpoint directory (``Geneformer-V2-104M``
    from https://huggingface.co/ctheodoris/Geneformer). Set
    ``model_directory`` to the path of that checkpoint.

    Embedding extraction is CLS-token based, matching Geneformer's default
    ``emb_mode='cls'``. Perturbation is implemented by setting the gene's
    raw-count column to zero (knockout) or a saturating value (overexpress)
    before rank-tokenization — equivalent in effect to the official
    ``InSilicoPerturber(perturb_type='delete'|'activate')`` because rank
    tokenization rank-orders genes by normalised expression, so a zero'd
    gene drops out of the top-N token sequence and a saturated gene sits at
    rank 1.

    The backend is stateful: a first ``embed()`` call lazy-loads the
    Geneformer tokenizer dictionaries + HuggingFace BERT model into RAM.
    Subsequent calls reuse the in-memory model.
    """

    # Geneformer V2 104M hidden size. Overwritten at load-time if the
    # checkpoint reports something different.
    embedding_dim = 768

    def __init__(
        self,
        model_directory: Optional[str] = None,
        device: str = "cpu",
        forward_batch_size: int = 8,
        max_input_len: int = 4096,
        emb_mode: str = "cls",
        cache_dir: Optional[str] = None,
    ) -> None:
        if emb_mode not in ("cls", "mean"):
            raise ValueError("emb_mode must be 'cls' or 'mean'")
        self.model_directory = model_directory
        self.device = device
        self.forward_batch_size = int(forward_batch_size)
        self.max_input_len = int(max_input_len)
        self.emb_mode = emb_mode
        # Optional content-hashed on-disk cache for ``embed()`` results.
        # The key covers (model_directory basename, emb_mode, max_input_len,
        # expression bytes), so zeroed-gene perturbation inputs produce a
        # distinct key from the baseline automatically. Reruns with the
        # same cohort hit the cache and skip the 40-minute forward pass.
        self.cache_dir: Optional[Path] = Path(cache_dir).expanduser() if cache_dir else None
        self._model: Any = None
        self._torch: Any = None
        self._gene_median: Optional[Dict[str, float]] = None
        self._token_dict: Optional[Dict[str, int]] = None
        self._gene_mapping: Optional[Dict[str, str]] = None
        self._cls_token_id: Optional[int] = None
        self._eos_token_id: Optional[int] = None
        self._pad_token_id: Optional[int] = None

    # --------------------------------------------------------------- lazy ---
    def _lazy_load(self) -> None:
        if self._model is not None:
            return
        try:
            import pickle  # stdlib
            from pathlib import Path

            import torch
            from transformers import BertForMaskedLM, BertModel
        except ImportError as exc:  # pragma: no cover - optional dep
            raise ImportError(
                "OfficialBackend requires torch + transformers + geneformer. "
                "Install with: pip install pathway-subtyping[perturb]"
            ) from exc

        import geneformer  # type: ignore

        self._torch = torch
        geneformer_pkg = Path(geneformer.__file__).parent

        with open(geneformer_pkg / "gene_median_dictionary_gc104M.pkl", "rb") as f:
            self._gene_median = pickle.load(f)
        with open(geneformer_pkg / "token_dictionary_gc104M.pkl", "rb") as f:
            self._token_dict = pickle.load(f)
        with open(geneformer_pkg / "ensembl_mapping_dict_gc104M.pkl", "rb") as f:
            self._gene_mapping = pickle.load(f)

        self._cls_token_id = self._token_dict["<cls>"]
        self._eos_token_id = self._token_dict["<eos>"]
        self._pad_token_id = self._token_dict["<pad>"]

        model_dir = self.model_directory
        if model_dir is None:
            raise ValueError(
                "OfficialBackend requires model_directory; set it to the "
                "path of a Geneformer V2 checkpoint (e.g., the "
                "'Geneformer-V2-104M' directory from the HuggingFace repo)."
            )
        # Load the BertModel backbone (dropping the LM head — we only need
        # hidden states). The Geneformer checkpoint is BertForMaskedLM; we
        # instantiate BertForMaskedLM and then use .bert for embeddings.
        model = BertForMaskedLM.from_pretrained(str(model_dir))
        model.eval()
        model.to(self.device)
        self._model = model.bert  # pure BertModel for hidden states
        # Reconcile embedding_dim with the actual checkpoint
        self.embedding_dim = int(model.config.hidden_size)

    # ----------------------------------------------------- preprocessing ---
    def _symbols_to_ensembl(
        self,
        symbols: List[str],
    ) -> Tuple[List[str], List[int]]:
        """Map a list of gene symbols to Ensembl IDs, filtering unmapped.

        Returns (ensembl_ids, keep_indices) where keep_indices are the
        positions in ``symbols`` that mapped and appear in the token
        dictionary.
        """
        assert self._gene_mapping is not None
        assert self._token_dict is not None
        keep: List[int] = []
        ens: List[str] = []
        for i, sym in enumerate(symbols):
            eid = self._gene_mapping.get(sym)
            if eid is None:
                continue
            if eid not in self._token_dict:
                continue
            keep.append(i)
            ens.append(eid)
        return ens, keep

    def _rank_tokenize_one(
        self,
        raw_counts: np.ndarray,
        ensembl_ids: List[str],
        target_sum: float = 10_000.0,
    ) -> List[int]:
        """Rank-based tokenization mirroring Geneformer's tokenizer.

        Given raw counts for one cell and the corresponding Ensembl IDs,
        compute per-gene median-normalised rank and return the token IDs
        of the top-``max_input_len - 2`` genes ordered by rank, sandwiched
        between <cls> ... <eos>.
        """
        assert self._gene_median is not None
        assert self._token_dict is not None
        assert self._cls_token_id is not None and self._eos_token_id is not None

        raw = np.asarray(raw_counts, dtype=float)
        n = raw.sum()
        if n <= 0:
            return [self._cls_token_id, self._eos_token_id]
        norm = raw / n * target_sum
        medians = np.asarray(
            [self._gene_median.get(eid, np.nan) for eid in ensembl_ids],
        )
        valid = ~np.isnan(medians) & (medians > 0) & (norm > 0)
        if not valid.any():
            return [self._cls_token_id, self._eos_token_id]
        scaled = np.full_like(norm, 0.0)
        scaled[valid] = norm[valid] / medians[valid]
        # Rank by decreasing scaled expression, keep top-K
        order = np.argsort(-scaled)
        # Drop zero-scaled entries entirely
        order = order[scaled[order] > 0]
        token_budget = self.max_input_len - 2  # leave space for cls/eos
        order = order[:token_budget]
        tokens = [self._token_dict[ensembl_ids[i]] for i in order]
        return [self._cls_token_id] + tokens + [self._eos_token_id]

    def _tokenize_df(
        self,
        expression: pd.DataFrame,
    ) -> Tuple[List[List[int]], List[int]]:
        """Convert a samples x genes DataFrame to per-sample token sequences.

        Returns (token_sequences, kept_sample_indices). Samples whose entire
        expression vector zeros out (empty token sequence) are dropped.
        """
        if not isinstance(expression, pd.DataFrame):
            raise TypeError("OfficialBackend requires a DataFrame input")
        self._lazy_load()
        symbols = list(expression.columns)
        ensembl_ids, keep_gene_idx = self._symbols_to_ensembl(symbols)
        if not ensembl_ids:
            raise RuntimeError(
                "no gene symbols mapped to Geneformer vocabulary; check "
                "that expression columns are HGNC symbols"
            )
        X = expression.to_numpy(dtype=float)[:, keep_gene_idx]

        token_seqs: List[List[int]] = []
        kept_samples: List[int] = []
        for i in range(X.shape[0]):
            seq = self._rank_tokenize_one(X[i], ensembl_ids)
            if len(seq) <= 2:  # only cls+eos — zero expression
                continue
            token_seqs.append(seq)
            kept_samples.append(i)
        return token_seqs, kept_samples

    # ------------------------------------------------------ forward pass ---
    def _forward(self, token_seqs: List[List[int]]) -> np.ndarray:
        """Batched forward pass; returns (n, embedding_dim) array."""
        torch = self._torch
        assert self._pad_token_id is not None
        out = np.zeros((len(token_seqs), self.embedding_dim), dtype=float)
        for start in range(0, len(token_seqs), self.forward_batch_size):
            batch = token_seqs[start : start + self.forward_batch_size]
            max_len = max(len(s) for s in batch)
            input_ids = np.full(
                (len(batch), max_len),
                self._pad_token_id,
                dtype=np.int64,
            )
            attn = np.zeros((len(batch), max_len), dtype=np.int64)
            for i, seq in enumerate(batch):
                input_ids[i, : len(seq)] = seq
                attn[i, : len(seq)] = 1
            input_ids_t = torch.as_tensor(input_ids, device=self.device)
            attn_t = torch.as_tensor(attn, device=self.device)
            with torch.no_grad():
                res = self._model(
                    input_ids=input_ids_t,
                    attention_mask=attn_t,
                )
            hidden = res.last_hidden_state.cpu().numpy()  # (b, L, H)
            if self.emb_mode == "cls":
                emb = hidden[:, 0, :]
            else:  # mean-pool over non-padding tokens
                mask = attn[:, :, None].astype(float)
                emb = (hidden * mask).sum(axis=1) / mask.sum(axis=1).clip(1e-9)
            out[start : start + len(batch)] = emb
        return out

    # ----------------------------------------------- content-hash cache ---
    def _backend_id(self) -> str:
        """Stable identifier that scopes cache keys to checkpoint + config."""
        checkpoint_tag = Path(self.model_directory).name if self.model_directory else "none"
        return (
            f"geneformer:official:{checkpoint_tag}" f":{self.emb_mode}:maxlen{self.max_input_len}"
        )

    def _cache_lookup(
        self,
        expression: pd.DataFrame,
    ) -> Optional[np.ndarray]:
        if self.cache_dir is None:
            return None
        from pathway_subtyping.embed.cache import cache_key_for

        key = cache_key_for(
            backend=self._backend_id(),
            expression=expression,
        )
        path = self.cache_dir / f"{key}.npy"
        if path.exists():
            logger.info(
                "[OfficialBackend] cache hit: %s… n=%d",
                key[:12],
                len(expression),
            )
            return np.load(path)
        return None

    def _cache_store(
        self,
        expression: pd.DataFrame,
        arr: np.ndarray,
    ) -> None:
        if self.cache_dir is None:
            return
        from pathway_subtyping.embed.cache import cache_key_for

        key = cache_key_for(
            backend=self._backend_id(),
            expression=expression,
        )
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        np.save(self.cache_dir / f"{key}.npy", arr)
        logger.info(
            "[OfficialBackend] cache stored: %s… n=%d d=%d",
            key[:12],
            arr.shape[0],
            arr.shape[1],
        )

    # ------------------------------------------------------ public API ---
    def embed(self, expression: pd.DataFrame) -> np.ndarray:
        cached = self._cache_lookup(expression)
        if cached is not None:
            return cached
        self._lazy_load()
        seqs, kept = self._tokenize_df(expression)
        full = np.zeros((len(expression), self.embedding_dim), dtype=float)
        if seqs:
            full[kept] = self._forward(seqs)
        self._cache_store(expression, full)
        return full

    def perturb(
        self,
        expression: pd.DataFrame,
        gene: str,
        mode: PerturbationMode,
    ) -> np.ndarray:
        if gene not in expression.columns:
            raise KeyError(
                f"gene {gene!r} not in expression columns ({len(expression.columns)} cols)"
            )
        perturbed = expression.copy()
        if mode == PerturbationMode.KNOCKOUT:
            perturbed[gene] = 0.0
        elif mode == PerturbationMode.OVEREXPRESS:
            # Saturate to 10x the per-cell maximum — guarantees rank 1.
            col_max = perturbed.max(axis=1)
            perturbed[gene] = col_max * 10.0
        else:
            raise ValueError(f"unknown mode: {mode}")
        return self.embed(perturbed)


# --------------------------------------------------------------------------- #
# Deterministic fallback
# --------------------------------------------------------------------------- #


class FallbackPerturber(GeneformerBackend):
    """Deterministic PCA-based substitute for Geneformer.

    Fits a PCA basis on the first expression matrix it sees. ``embed``
    projects into that basis; ``perturb(expr, gene, mode)`` simulates
    a perturbation by setting the gene's column to zero (knockout) or
    amplifying it (overexpress) before projecting.

    This is an API harness, not a biology-grade substitute. The
    directional behaviour on well-separated markers is correct —
    knocking out a gene that anchors a cluster does move cells toward
    the mean embedding — but absolute magnitudes are meaningless.
    """

    def __init__(
        self,
        embedding_dim: int = 64,
        overexpression_factor: float = 3.0,
        seed: Optional[int] = 0,
    ) -> None:
        if embedding_dim < 2:
            raise ValueError("embedding_dim must be >= 2")
        self.embedding_dim = int(embedding_dim)
        self.overexpression_factor = float(overexpression_factor)
        self.seed = seed
        self._mean: Optional[np.ndarray] = None
        self._components: Optional[np.ndarray] = None
        self._columns: Optional[List[str]] = None
        self._fitted: bool = False

    # --------------------------------------------------------------- fit ---
    def fit(self, expression: pd.DataFrame) -> "FallbackPerturber":
        if not isinstance(expression, pd.DataFrame):
            raise TypeError("FallbackPerturber requires a DataFrame input")
        self._columns = list(expression.columns)
        X = expression.to_numpy(dtype=float)
        self._mean = X.mean(axis=0)
        centered = X - self._mean
        k = min(self.embedding_dim, centered.shape[1], centered.shape[0])
        _, _, Vt = np.linalg.svd(centered, full_matrices=False)
        self._components = Vt[:k]
        if k < self.embedding_dim:
            rng = np.random.default_rng(self.seed)
            extra = rng.standard_normal((self.embedding_dim - k, centered.shape[1]))
            extra = np.linalg.qr(extra.T)[0].T
            self._components = np.vstack([self._components, extra])[: self.embedding_dim]
        self._fitted = True
        logger.info(
            "[FallbackPerturber] fitted: embedding_dim=%d n_genes=%d",
            self.embedding_dim,
            centered.shape[1],
        )
        return self

    # -------------------------------------------------------------- embed ---
    def embed(self, expression: pd.DataFrame) -> np.ndarray:
        if not self._fitted:
            self.fit(expression)
        assert self._mean is not None and self._components is not None
        X = expression.to_numpy(dtype=float)
        if X.shape[1] != self._components.shape[1]:
            raise ValueError(
                f"expression has {X.shape[1]} genes; backend was fitted on "
                f"{self._components.shape[1]}"
            )
        return (X - self._mean) @ self._components.T

    # ------------------------------------------------------------ perturb ---
    def perturb(
        self,
        expression: pd.DataFrame,
        gene: str,
        mode: PerturbationMode,
    ) -> np.ndarray:
        if not self._fitted:
            self.fit(expression)
        if gene not in expression.columns:
            raise KeyError(
                f"gene {gene!r} not in expression columns; got "
                f"{len(expression.columns)} columns"
            )
        perturbed = expression.copy()
        if mode == PerturbationMode.KNOCKOUT:
            perturbed[gene] = 0.0
        elif mode == PerturbationMode.OVEREXPRESS:
            perturbed[gene] = perturbed[gene] * self.overexpression_factor
        else:
            raise ValueError(f"unknown mode: {mode}")
        return self.embed(perturbed)


# --------------------------------------------------------------------------- #
# High-level wrapper
# --------------------------------------------------------------------------- #


@dataclass
class PerturbationResult:
    """One perturbation call's result.

    Attributes:
        gene: The gene that was perturbed.
        mode: Knockout or overexpress.
        baseline_embedding: Embedding of the unperturbed expression (n_cells, d).
        perturbed_embedding: Embedding after perturbation (n_cells, d).
        delta_embedding: ``perturbed_embedding - baseline_embedding``.
        l2_impact_per_cell: L2 norm of delta per cell.
    """

    gene: str
    mode: PerturbationMode
    baseline_embedding: np.ndarray
    perturbed_embedding: np.ndarray
    delta_embedding: np.ndarray
    l2_impact_per_cell: np.ndarray

    @property
    def mean_l2_impact(self) -> float:
        return float(self.l2_impact_per_cell.mean())

    def to_dict(self) -> Dict[str, Any]:
        return {
            "gene": self.gene,
            "mode": self.mode.value,
            "n_cells": int(self.baseline_embedding.shape[0]),
            "embedding_dim": int(self.baseline_embedding.shape[1]),
            "mean_l2_impact": self.mean_l2_impact,
        }


class GeneformerPerturber:
    """High-level perturbation interface.

    Usage::

        perturber = GeneformerPerturber()  # uses FallbackPerturber by default
        result = perturber.perturb(expression_df, gene="MYC", mode="knockout")
        # result.delta_embedding is (n_cells, embedding_dim)
    """

    def __init__(
        self,
        backend: Optional[GeneformerBackend] = None,
    ) -> None:
        self.backend = backend if backend is not None else FallbackPerturber()

    @property
    def embedding_dim(self) -> int:
        return self.backend.embedding_dim

    # -------------------------------------------------------------- embed ---
    def embed(self, expression: pd.DataFrame) -> np.ndarray:
        return self.backend.embed(expression)

    # ----------------------------------------------------------- perturb ---
    def perturb(
        self,
        expression: pd.DataFrame,
        gene: str,
        mode: PerturbationMode = PerturbationMode.KNOCKOUT,
    ) -> PerturbationResult:
        mode = PerturbationMode(mode) if isinstance(mode, str) else mode
        baseline = self.backend.embed(expression)
        perturbed = self.backend.perturb(expression, gene, mode)
        delta = perturbed - baseline
        l2 = np.linalg.norm(delta, axis=1)
        return PerturbationResult(
            gene=gene,
            mode=mode,
            baseline_embedding=baseline,
            perturbed_embedding=perturbed,
            delta_embedding=delta,
            l2_impact_per_cell=l2,
        )
