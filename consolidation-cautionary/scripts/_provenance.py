"""Shared provenance recording for the reproduction packages.

Why this exists
---------------
Every package that fetches from a live service (cBioPortal, NCBI E-utilities) was
previously depositing a result JSON that named, at best, a study id. That is not enough
to tie a number to a data snapshot. cBioPortal is not a versioned API: a study can gain
samples, drop samples, or change its molecular profiles between one fetch and the next,
and nothing in the deposited artifact would reveal that. A reviewer who re-fetched and
got different numbers had no way to tell an upstream change from an analysis error.

The blocks below are recorded into every network-derived result JSON so that question is
always answerable:

* ``env_provenance()``   — the library versions that can move a number. torch is included
  because Result 3's DL baselines are not deterministic across torch releases, and the
  manuscript quotes a specific version.
* ``fetch_provenance()`` — endpoint, study, UTC fetch date, post-filter shape, and a
  SHA-256 of the assembled matrix. The hash is the load-bearing part: two fetches that
  agree on it are working from identical data, whatever the date says.

Deliberately NOT recorded: hostname, username, absolute paths. These artifacts are
published, and machine identity is not provenance.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import importlib
import platform


def cache_write(df, path: str) -> str:
    """Deposit an assembled input matrix so the analysis becomes offline-reproducible.

    Written with ``mtime=0`` because gzip stamps the current time into its header, which
    would make the identical matrix hash differently on every run and produce spurious
    MANIFEST mismatches for a reviewer. Returns the SHA-256 of the contents.
    """
    import os

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    df.to_csv(path, compression={"method": "gzip", "mtime": 0})
    return sha256_frame(df)


def cache_read(path: str):
    """Read a deposited input matrix written by :func:`cache_write`."""
    import pandas as pd

    return pd.read_csv(path, index_col=0)


def _ver(mod: str):
    """Return an installed module's version, or None if it is absent."""
    try:
        return getattr(importlib.import_module(mod), "__version__", "unknown")
    except Exception:
        return None


def env_provenance() -> dict:
    """Library versions that can move a deposited number."""
    env = {
        "python": platform.python_version(),
        "pathway_subtyping": _ver("pathway_subtyping"),
        "numpy": _ver("numpy"),
        "pandas": _ver("pandas"),
        "scipy": _ver("scipy"),
        "sklearn": _ver("sklearn"),
        "statsmodels": _ver("statsmodels"),
        # torch is optional: only the DL baselines need it. Recorded as None when absent
        # so a reviewer can tell "not used" apart from "not recorded".
        "torch": _ver("torch"),
    }
    env["note"] = (
        "scikit-learn is the reproduction-sensitive pin (GMM/BIC determinism). "
        "torch is null unless the DL baselines ran; those are NOT deterministic "
        "across torch releases, so compare this field before reporting a mismatch."
    )
    return env


def sha256_frame(df) -> str:
    """Stable SHA-256 of a DataFrame's contents, independent of file formatting.

    Hashes the CSV serialisation rather than the in-memory buffer so that dtype
    promotion (int64 vs float64 on a re-fetch) does not change the digest for
    numerically identical data.
    """
    return hashlib.sha256(df.to_csv(index=True).encode("utf-8")).hexdigest()


def fetch_provenance(
    endpoint: str,
    study: str,
    matrix=None,
    n_samples: int | None = None,
    n_features: int | None = None,
    extra: dict | None = None,
) -> dict:
    """Record what was fetched, from where, when, and what it hashed to.

    Args:
        endpoint: Base URL of the service (no credentials, no query string).
        study:    Study / accession identifier as the service names it.
        matrix:   Optional assembled DataFrame; hashed and shaped if given.
        n_samples, n_features: Explicit counts when no matrix is passed.
        extra:    Package-specific fields (filters applied, profile ids, ...).
    """
    prov = {
        "endpoint": endpoint,
        "study": study,
        # Date, not datetime: the useful granularity for detecting upstream drift, and
        # it keeps the artifact stable across same-day re-runs.
        "fetch_date_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d"),
        "service_is_versioned": False,
        "drift_note": (
            "This endpoint is NOT versioned. If your re-run disagrees, compare "
            "matrix_sha256 first: a different hash means the upstream data changed, "
            "not that the analysis did."
        ),
    }
    if matrix is not None:
        prov["n_samples"] = int(matrix.shape[0])
        prov["n_features"] = int(matrix.shape[1])
        prov["matrix_sha256"] = sha256_frame(matrix)
    if n_samples is not None:
        prov["n_samples"] = int(n_samples)
    if n_features is not None:
        prov["n_features"] = int(n_features)
    if extra:
        prov.update(extra)
    return prov
