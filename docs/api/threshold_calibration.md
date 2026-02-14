# Threshold Calibration API Reference

The `threshold_calibration` module provides data-driven validation threshold calibration, replacing hard-coded thresholds with values that adjust for sample size and number of clusters.

## Overview

Hard-coded thresholds (0.15 for null ARI, 0.8 for stability) don't account for:
- **Sample size**: Small samples produce noisier ARI distributions
- **Cluster count**: More clusters inflate chance ARI

This module uses pre-computed lookup tables derived from empirical simulations to provide appropriate thresholds for any (n_samples, n_clusters) configuration.

## Classes

### CalibratedThresholds

Calibrated validation thresholds for a specific dataset configuration.

```python
from pathway_subtyping import CalibratedThresholds
```

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `null_ari_threshold` | `float` | Maximum ARI under null hypothesis |
| `stability_threshold` | `float` | Minimum ARI for stability test |
| `n_samples` | `int` | Number of samples calibrated for |
| `n_clusters` | `int` | Number of clusters calibrated for |
| `n_pathways` | `int` | Number of pathways used |
| `alpha` | `float` | Significance level |
| `calibration_method` | `str` | Method used: "lookup", "interpolated", "simulated", "default" |
| `interpolated` | `bool` | Whether values were interpolated |

#### Methods

##### `to_dict() -> Dict[str, Any]`

Convert to dictionary for JSON serialization.

##### `format_report() -> str`

Format as human-readable report string.

```python
ct = calibrate_thresholds(n_samples=100, n_clusters=3)
print(ct.format_report())
```

##### `get_citations() -> List[str]`

Return relevant citations for the calibration methodology.

---

### CalibrationSimulationResult

Results from a calibration simulation run.

```python
from pathway_subtyping import CalibrationSimulationResult
```

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `null_ari_distribution` | `np.ndarray` | Distribution of null ARI values |
| `stability_distribution` | `np.ndarray` | Distribution of stability ARI values |
| `null_ari_percentile` | `float` | Percentile value for null threshold |
| `stability_percentile` | `float` | Percentile value for stability threshold |
| `n_simulations` | `int` | Number of simulations run |

#### Methods

##### `to_dict() -> Dict[str, Any]`

Convert to dictionary for JSON serialization.

---

## Functions

### `calibrate_thresholds(...) -> CalibratedThresholds`

Calibrate validation thresholds for a given dataset configuration.

```python
from pathway_subtyping import calibrate_thresholds

ct = calibrate_thresholds(
    n_samples=150,
    n_clusters=4,
    n_pathways=15,
    alpha=0.05,
    method="auto",
    n_simulations=200,
    seed=42,
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `n_samples` | `int` | required | Number of samples in dataset |
| `n_clusters` | `int` | required | Number of clusters |
| `n_pathways` | `int` | `15` | Number of pathways |
| `alpha` | `float` | `0.05` | Significance level |
| `method` | `str` | `"auto"` | Calibration method: "auto", "lookup", "simulate" |
| `n_simulations` | `int` | `200` | Simulations for "simulate" method |
| `seed` | `Optional[int]` | `None` | Random seed for reproducibility |

**Methods:**
- `"auto"`: Try lookup first, interpolate if needed, simulate if out of range
- `"lookup"`: Lookup/interpolate only; raises `ValueError` if out of range
- `"simulate"`: Always run fresh simulations (slower but exact)

**Returns:** `CalibratedThresholds`

---

### `get_default_thresholds() -> CalibratedThresholds`

Return default (legacy) thresholds for backward compatibility.

```python
from pathway_subtyping import get_default_thresholds

defaults = get_default_thresholds()
# defaults.null_ari_threshold == 0.15
# defaults.stability_threshold == 0.8
```

---

### `generate_calibration_table(...) -> Dict[str, Dict]`

Generate calibration lookup tables from simulations. Long-running; used to regenerate the pre-computed tables.

```python
from pathway_subtyping import generate_calibration_table

tables = generate_calibration_table(
    sample_sizes=[50, 100, 200],
    cluster_counts=[2, 3, 4],
    n_simulations=500,
    seed=42,
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sample_sizes` | `Optional[List[int]]` | grid default | Sample sizes to simulate |
| `cluster_counts` | `Optional[List[int]]` | grid default | Cluster counts to simulate |
| `n_simulations` | `int` | `500` | Simulations per grid point |
| `seed` | `Optional[int]` | `None` | Random seed |

**Returns:** Dict with `"null_ari"` and `"stability"` sub-dicts mapping `(n_samples, n_clusters)` tuples to threshold values.

---

## Pre-Computed Lookup Tables

The module includes pre-computed tables covering:

- **Sample sizes**: 30, 50, 75, 100, 150, 200, 300, 500
- **Cluster counts**: 2, 3, 4, 5, 6, 7, 8
- **Total grid points**: 56

**Properties:**
- Null ARI thresholds *decrease* with more samples (tighter null distribution)
- Null ARI thresholds *increase* with more clusters (more chance agreement)
- Stability thresholds *increase* with more samples (more stable clustering)
- Stability thresholds *decrease* with more clusters (harder to maintain stability)

Tables were generated with `n_simulations=500` and `seed=42`.

---

## Example: Pipeline Integration

The pipeline automatically uses threshold calibration when thresholds are `null`:

```yaml
validation:
  run_gates: true
  calibrate: true
  stability_threshold: null   # Auto-calibrated
  null_ari_max: null           # Auto-calibrated
```

Or programmatically:

```python
from pathway_subtyping import calibrate_thresholds
from pathway_subtyping.validation import ValidationGates

# Calibrate for your data
ct = calibrate_thresholds(n_samples=200, n_clusters=4)

# Use calibrated thresholds
validator = ValidationGates(
    stability_threshold=ct.stability_threshold,
    null_ari_max=ct.null_ari_threshold,
)
```

---

## Regenerating Lookup Tables

```bash
python scripts/generate_calibration_table.py --n_simulations 500 --seed 42
```

This takes ~10-30 minutes and prints Python dict literals for copy-paste into the module.
