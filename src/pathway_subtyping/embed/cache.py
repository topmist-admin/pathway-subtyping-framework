"""
Content-hashed embedding cache.

The cache key covers the backend identifier (which includes model
checkpoint version), the tokenizer / config version (if any), and the
bytes of the input expression matrix. Changing any of those
automatically invalidates the cached entry — no manual version bumps
needed.

Usage::

    cache = EmbeddingCache(cache_dir="~/.cache/pathway-subtyping/embed")
    key = cache_key_for(backend="scgpt:official/v1.2.0", expression=expr)
    if cache.has(key):
        result = cache.get(key)
    else:
        result = embedder.embed(expr)
        cache.put(key, result)

The cache is safe across package upgrades: if a new PSF release changes
the backend identifier or tokenizer version, cached entries produced by
the old version will no longer match the key and a fresh inference pass
is triggered.

Research use only. Not for clinical decision-making.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

from .base import EmbeddingResult

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# Hashing
# --------------------------------------------------------------------------- #

def _hash_expression(expression: pd.DataFrame) -> str:
    """Deterministic hash over the expression content.

    Covers column order, row order, and numeric values. Uses float64
    byte representation so floats compare exactly across platforms that
    honour IEEE-754.
    """
    h = hashlib.sha256()
    h.update(str(list(expression.columns)).encode("utf-8"))
    h.update(b"|")
    h.update(str(list(expression.index)).encode("utf-8"))
    h.update(b"|")
    arr = np.ascontiguousarray(
        expression.to_numpy(dtype=float), dtype=np.float64
    )
    h.update(arr.tobytes())
    return h.hexdigest()


def cache_key_for(
    backend: str,
    expression: pd.DataFrame,
    extra: Optional[Dict[str, Any]] = None,
) -> str:
    """Build a single cache key from backend + expression + optional extras.

    Args:
        backend: Short backend identifier, e.g. ``'scgpt:official/v1.2.0'``.
            Changing checkpoints → different backend string → cache miss.
        expression: Input DataFrame whose content must be part of the key.
        extra: Optional dict of extra key material (e.g. tokenizer args).
    """
    h = hashlib.sha256()
    h.update(backend.encode("utf-8"))
    h.update(b"|")
    h.update(_hash_expression(expression).encode("utf-8"))
    if extra is not None:
        h.update(b"|")
        h.update(
            json.dumps(extra, sort_keys=True, default=str).encode("utf-8")
        )
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# On-disk cache
# --------------------------------------------------------------------------- #

class EmbeddingCache:
    """Simple filesystem cache for :class:`EmbeddingResult` objects.

    Stores each entry as two files: a numpy ``.npy`` for the embedding
    matrix and a JSON sidecar with the cell index, backend, and meta.
    """

    def __init__(self, cache_dir: str | Path):
        self.cache_dir = Path(cache_dir).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------- helpers ---
    def _paths(self, key: str) -> tuple[Path, Path]:
        return (
            self.cache_dir / f"{key}.npy",
            self.cache_dir / f"{key}.json",
        )

    def has(self, key: str) -> bool:
        npy, meta = self._paths(key)
        return npy.exists() and meta.exists()

    # --------------------------------------------------------------- get ---
    def get(self, key: str) -> EmbeddingResult:
        if not self.has(key):
            raise KeyError(f"no cached embedding for key {key[:12]}…")
        npy, meta_path = self._paths(key)
        arr = np.load(npy)
        with meta_path.open("r") as fh:
            meta_json = json.load(fh)
        logger.info("[EmbeddingCache] hit: %s…", key[:12])
        return EmbeddingResult(
            embeddings=arr,
            cell_index=pd.Index(meta_json["cell_index"], name=meta_json.get("index_name")),
            backend=meta_json["backend"],
            meta=meta_json.get("meta", {}),
        )

    # --------------------------------------------------------------- put ---
    def put(self, key: str, result: EmbeddingResult) -> None:
        npy, meta_path = self._paths(key)
        np.save(npy, result.embeddings)
        payload = {
            "backend": result.backend,
            "cell_index": list(result.cell_index),
            "index_name": result.cell_index.name,
            "meta": result.meta,
        }
        with meta_path.open("w") as fh:
            json.dump(payload, fh, default=str)
        logger.info(
            "[EmbeddingCache] stored: %s… (n=%d, d=%d)",
            key[:12], result.n_cells, result.embedding_dim,
        )

    def delete(self, key: str) -> bool:
        """Remove a cached entry. Returns True iff anything was deleted."""
        npy, meta_path = self._paths(key)
        deleted = False
        if npy.exists():
            npy.unlink()
            deleted = True
        if meta_path.exists():
            meta_path.unlink()
            deleted = True
        return deleted

    def clear(self) -> int:
        """Remove every cached entry under ``cache_dir``. Returns count removed."""
        n = 0
        for p in self.cache_dir.iterdir():
            if p.suffix in {".npy", ".json"}:
                p.unlink()
                n += 1
        return n
