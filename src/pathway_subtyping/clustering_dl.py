"""Deep-learning clustering baselines (DEC, VAE-GMM).

Added to answer Scientific Reports reviewer point R2.2: the manuscript compared
pathway-level clustering against PCA/spectral/Leiden/k-means/NMF but not modern
deep-learning clustering. This module adds two standard DL clusterers so they can
be (a) benchmarked as methods and (b) — the load-bearing point — run through the
same validation gates, demonstrating that the gates are *clusterer-agnostic*: a
DL method's clusters need discreteness validation exactly like a GMM's.

Both functions match the comparison-harness convention used by the sklearn
baselines: ``run_x(X, k, ...) -> (labels: np.ndarray, silhouette: float)``.

torch is imported lazily so importing this module never fails; the functions
raise a clear error if torch is unavailable.

References
----------
- DEC: Xie, Girshick & Farhadi, "Unsupervised Deep Embedding for Clustering
  Analysis", ICML 2016. arXiv:1511.06335
- VaDE / VAE-GMM: Jiang et al., "Variational Deep Embedding: An Unsupervised and
  Generative Approach to Clustering", IJCAI 2017. arXiv:1611.05148
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler


def _require_torch():
    try:
        import torch  # noqa: F401
        return torch
    except ImportError as e:  # pragma: no cover - environment-dependent
        raise ImportError(
            "clustering_dl requires torch. Install with `pip install torch` "
            "(already present in the PSF dev environment)."
        ) from e


def _latent_dim(n: int, p: int) -> int:
    """Small latent — deliberately narrow at the small n this framework targets."""
    return int(max(2, min(10, p // 5, n // 5)))


def _safe_silhouette(X: np.ndarray, labels: np.ndarray) -> float:
    if len(set(labels)) < 2:
        return 0.0
    try:
        return float(silhouette_score(X, labels))
    except Exception:
        return 0.0


def _build_mlp(torch, dims):
    """Sequential MLP with ReLU between layers (no activation on the last)."""
    nn = torch.nn
    layers = []
    for i in range(len(dims) - 1):
        layers.append(nn.Linear(dims[i], dims[i + 1]))
        if i < len(dims) - 2:
            layers.append(nn.ReLU())
    return nn.Sequential(*layers)


def run_dec(
    X,
    k: int,
    seed: int = 42,
    latent_dim: Optional[int] = None,
    pretrain_epochs: int = 150,
    cluster_epochs: int = 60,
    lr: float = 1e-3,
) -> Tuple[np.ndarray, float]:
    """Deep Embedded Clustering (Xie et al. 2016).

    Pretrain an autoencoder, initialise cluster centroids with k-means on the
    latent code, then jointly refine the embedding and assignments by minimising
    KL(target || soft-assignment). Returns (labels, silhouette-in-input-space).
    """
    torch = _require_torch()
    torch.manual_seed(seed)
    Xz = StandardScaler().fit_transform(np.asarray(X, dtype=np.float64))
    n, p = Xz.shape
    d = latent_dim or _latent_dim(n, p)
    Xt = torch.tensor(Xz, dtype=torch.float32)

    enc = _build_mlp(torch, [p, max(32, 2 * d), d])
    dec = _build_mlp(torch, [d, max(32, 2 * d), p])
    params = list(enc.parameters()) + list(dec.parameters())
    opt = torch.optim.Adam(params, lr=lr)
    mse = torch.nn.MSELoss()

    # 1) pretrain autoencoder (reconstruction)
    for _ in range(pretrain_epochs):
        opt.zero_grad()
        z = enc(Xt)
        loss = mse(dec(z), Xt)
        loss.backward()
        opt.step()

    # 2) init centroids via k-means on the latent code
    with torch.no_grad():
        Z0 = enc(Xt).numpy()
    km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(Z0)
    centroids = torch.tensor(km.cluster_centers_, dtype=torch.float32, requires_grad=True)

    opt2 = torch.optim.Adam(list(enc.parameters()) + [centroids], lr=lr)

    def soft_assign(z):
        # Student's t-distribution kernel (alpha=1), Xie et al. eq. 1
        dist2 = torch.cdist(z, centroids) ** 2
        q = 1.0 / (1.0 + dist2)
        return q / q.sum(1, keepdim=True)

    # 3) KL-sharpening refinement
    for _ in range(cluster_epochs):
        opt2.zero_grad()
        z = enc(Xt)
        q = soft_assign(z)
        # target distribution p (sharpened, normalised by cluster frequency)
        weight = q ** 2 / q.sum(0)
        p = (weight.t() / weight.sum(1)).t().detach()
        kl = (p * (torch.log(p + 1e-9) - torch.log(q + 1e-9))).sum(1).mean()
        kl.backward()
        opt2.step()

    with torch.no_grad():
        labels = soft_assign(enc(Xt)).argmax(1).numpy()
    return labels, _safe_silhouette(Xz, labels)


def run_vae_gmm(
    X,
    k: int,
    seed: int = 42,
    latent_dim: Optional[int] = None,
    epochs: int = 200,
    lr: float = 1e-3,
    beta: float = 1.0,
) -> Tuple[np.ndarray, float]:
    """VAE + GMM clustering (VaDE-style, Jiang et al. 2017).

    Train a variational autoencoder, then fit a Gaussian mixture on the latent
    means and assign. A pragmatic, stable stand-in for full VaDE that keeps the
    "deep generative embedding + mixture clustering" character. Returns
    (labels, silhouette-in-input-space).
    """
    torch = _require_torch()
    nn = torch.nn
    torch.manual_seed(seed)
    Xz = StandardScaler().fit_transform(np.asarray(X, dtype=np.float64))
    n, p = Xz.shape
    d = latent_dim or _latent_dim(n, p)
    Xt = torch.tensor(Xz, dtype=torch.float32)
    h = max(32, 2 * d)

    class VAE(nn.Module):
        def __init__(self):
            super().__init__()
            self.enc = nn.Sequential(nn.Linear(p, h), nn.ReLU())
            self.mu = nn.Linear(h, d)
            self.logvar = nn.Linear(h, d)
            self.dec = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Linear(h, p))

        def forward(self, x):
            he = self.enc(x)
            mu, logvar = self.mu(he), self.logvar(he)
            std = torch.exp(0.5 * logvar)
            z = mu + std * torch.randn_like(std)
            return self.dec(z), mu, logvar

    model = VAE()
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    for _ in range(epochs):
        opt.zero_grad()
        recon, mu, logvar = model(Xt)
        recon_loss = ((recon - Xt) ** 2).sum(1).mean()
        kld = -0.5 * (1 + logvar - mu ** 2 - logvar.exp()).sum(1).mean()
        (recon_loss + beta * kld).backward()
        opt.step()

    with torch.no_grad():
        Z = model.mu(model.enc(Xt)).numpy()
    gmm = GaussianMixture(n_components=k, covariance_type="full", n_init=10,
                          random_state=seed, reg_covar=1e-6).fit(Z)
    labels = gmm.predict(Z)
    return labels, _safe_silhouette(Xz, labels)


# Registry so harnesses can iterate uniformly.
DL_CLUSTERERS = {
    "dec": run_dec,
    "vae_gmm": run_vae_gmm,
}
