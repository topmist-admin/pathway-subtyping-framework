# Deep-Learning Clustering Baselines API

> **Module**: `pathway_subtyping.clustering_dl`
> **Added in**: v0.8.0

Two standard deep-clustering baselines — DEC and VAE-GMM (VaDE-style) — so that
deep-learning clusterers can be (a) benchmarked as methods against the pathway-GMM
pipeline and (b) run through the same validation gates, demonstrating that the gates
are *clusterer-agnostic*: a DL method's clusters need discreteness validation exactly
like a GMM's.

`torch` is imported **lazily**, so importing this module never fails. The functions
raise a clear `ImportError` if torch is unavailable.

---

## Quick Example

```python
from pathway_subtyping.clustering_dl import run_dec, run_vae_gmm

# X is a samples x features matrix (e.g. pathway_scores.values)
labels, sil = run_dec(X, k=3, seed=42)
print(labels, sil)

labels_v, sil_v = run_vae_gmm(X, k=3, seed=42)
```

Both functions follow the comparison-harness convention used by the sklearn
baselines:

```python
run_x(X, k, ...) -> (labels: np.ndarray, silhouette: float)
```

The returned silhouette is computed **in input space** (on the standardized input
matrix), not in the learned latent space, so it is directly comparable to the
silhouette reported by the non-DL baselines.

---

## Functions

### `run_dec()`

```python
def run_dec(
    X,
    k: int,
    seed: int = 42,
    latent_dim: Optional[int] = None,
    pretrain_epochs: int = 150,
    cluster_epochs: int = 60,
    lr: float = 1e-3,
) -> Tuple[np.ndarray, float]
```

Deep Embedded Clustering (Xie et al. 2016).

Pretrain an autoencoder, initialise cluster centroids with k-means on the latent
code, then jointly refine the embedding and assignments by minimising
KL(target || soft-assignment).

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `X` | array-like | required | Samples x features matrix (standardized internally) |
| `k` | int | required | Number of clusters |
| `seed` | int | 42 | Seed for `torch.manual_seed` and the k-means / centroid init |
| `latent_dim` | int or None | None | Latent dimension; when None, uses the internal small-latent rule |
| `pretrain_epochs` | int | 150 | Full-batch autoencoder reconstruction epochs |
| `cluster_epochs` | int | 60 | KL-sharpening refinement epochs |
| `lr` | float | 1e-3 | Adam learning rate (both phases) |

**Returns:** `(labels, silhouette)` — `labels` is the argmax of the soft assignment;
`silhouette` is the input-space silhouette (0.0 if fewer than 2 distinct labels or if
the silhouette computation fails).

**Raises:** `ImportError` if `torch` is not installed.

---

### `run_vae_gmm()`

```python
def run_vae_gmm(
    X,
    k: int,
    seed: int = 42,
    latent_dim: Optional[int] = None,
    epochs: int = 200,
    lr: float = 1e-3,
    beta: float = 1.0,
) -> Tuple[np.ndarray, float]
```

VAE + GMM clustering (VaDE-style, Jiang et al. 2017).

Train a variational autoencoder, then fit a Gaussian mixture on the latent means and
assign. A pragmatic, stable stand-in for full VaDE that keeps the "deep generative
embedding + mixture clustering" character.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `X` | array-like | required | Samples x features matrix (standardized internally) |
| `k` | int | required | Number of mixture components |
| `seed` | int | 42 | Seed for `torch.manual_seed` and `GaussianMixture(random_state=...)` |
| `latent_dim` | int or None | None | Latent dimension; when None, uses the internal small-latent rule |
| `epochs` | int | 200 | Full-batch VAE training epochs |
| `lr` | float | 1e-3 | Adam learning rate |
| `beta` | float | 1.0 | Weight on the KL divergence term of the ELBO |

**Returns:** `(labels, silhouette)` — `labels` from `GaussianMixture.predict` on the
latent means; `silhouette` is the input-space silhouette (0.0 if fewer than 2 distinct
labels or if the silhouette computation fails).

**Raises:** `ImportError` if `torch` is not installed.

---

## Registry

### `DL_CLUSTERERS`

```python
DL_CLUSTERERS = {
    "dec": run_dec,
    "vae_gmm": run_vae_gmm,
}
```

`Dict[str, Callable]` so harnesses can iterate uniformly over the DL baselines.

---

## Implementation Details

Both functions share the same preprocessing and sizing rules.

| Step | DEC | VAE-GMM |
|------|-----|---------|
| Preprocessing | `StandardScaler` on `float64` input | `StandardScaler` on `float64` input |
| Latent dim (when `latent_dim=None`) | `max(2, min(10, p // 5, n // 5))` | `max(2, min(10, p // 5, n // 5))` |
| Hidden width | `max(32, 2 * d)` | `max(32, 2 * d)` |
| Architecture | MLP encoder `p → max(32, 2d) → d`, decoder `d → max(32, 2d) → p`, ReLU between layers, no activation on the last layer | Encoder `Linear(p, h) + ReLU` with separate `mu` and `logvar` heads; decoder `Linear(d, h) + ReLU + Linear(h, p)` |
| Optimiser | Adam, full batch | Adam, full batch |
| Loss | Phase 1: MSE reconstruction. Phase 2: KL(p ‖ q) with a Student's t soft assignment (alpha = 1, Xie et al. eq. 1) and a frequency-normalised target `p` | Sum-over-features squared reconstruction error (mean over samples) + `beta` × KL divergence |
| Cluster assignment | Centroids initialised by `KMeans(n_clusters=k, n_init=10, random_state=seed)` on the latent code, then refined jointly with the encoder | `GaussianMixture(n_components=k, covariance_type="full", n_init=10, random_state=seed, reg_covar=1e-6)` on the latent means |

The latent dimension rule is deliberately narrow — the framework targets small *n*,
where a wide latent overfits.

---

## Optional Dependency

`torch` is not a core dependency. Both functions call an internal guard that raises:

```
ImportError: clustering_dl requires torch. Install with `pip install torch`
(already present in the PSF dev environment).
```

Because the import is lazy, `import pathway_subtyping.clustering_dl` itself always
succeeds — harnesses can import the module and skip the baselines gracefully when
torch is absent.

---

## References

- **DEC**: Xie, Girshick & Farhadi, "Unsupervised Deep Embedding for Clustering
  Analysis", ICML 2016. arXiv:1511.06335
- **VaDE / VAE-GMM**: Jiang et al., "Variational Deep Embedding: An Unsupervised and
  Generative Approach to Clustering", IJCAI 2017. arXiv:1611.05148

---

## See Also

- [Clustering](clustering.md) — the sklearn-based clustering algorithms and model selection
- [Discreteness Gate](discreteness.md) — clusterer-agnostic discreteness validation of any partition, including DL-produced ones
- [Benchmark](benchmark.md) — method comparison
