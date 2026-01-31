# Reproducibility Guide

> **RESEARCH USE ONLY** — This framework is for research purposes only. Not for clinical decision-making. See [DISCLAIMER.md](../DISCLAIMER.md).

This document describes the deterministic seed strategy and reproducibility guarantees for the Pathway Subtyping Framework.

## Seed Strategy

### Overview

All stochastic operations in the framework are controlled by a single master seed, ensuring fully reproducible results across runs.

### Configuration

Set the master seed in your config file:

```yaml
pipeline:
  seed: 42  # Any integer; null = random (non-reproducible)
```

**Important**: Always use an explicit seed for reproducible research.

### Seed Propagation

The master seed propagates to all randomized components:

| Component | Seeding Method | Notes |
|-----------|----------------|-------|
| Python random | `random.seed(seed)` | Standard library |
| NumPy | `np.random.seed(seed)` | Array operations |
| Scikit-learn | `random_state=seed` | GMM clustering |
| Bootstrap sampling | `np.random.RandomState(seed)` | Stability tests |

### Implementation

The framework uses a centralized seed utility:

```python
from pathway_subtyping.utils import set_global_seed, get_rng

# Set seed at pipeline entry point
set_global_seed(42)

# Get module-specific RNG for isolated reproducibility
rng = get_rng(42, module_name='clustering')
samples = rng.choice(data, size=100)
```

### Seed Utility Module

Location: `src/pathway_subtyping/utils/seed.py`

```python
"""Reproducibility utilities for deterministic execution."""

import random
import numpy as np
from typing import Optional

def set_global_seed(seed: int) -> None:
    """Set seed for all random number generators.

    Args:
        seed: Integer seed value
    """
    random.seed(seed)
    np.random.seed(seed)

def get_rng(seed: Optional[int], module_name: str = "") -> np.random.RandomState:
    """Get a module-specific random state for isolated reproducibility.

    Args:
        seed: Base seed (None for random)
        module_name: Module identifier for offset calculation

    Returns:
        NumPy RandomState instance
    """
    if seed is None:
        return np.random.RandomState()

    # Create deterministic offset from module name
    offset = sum(ord(c) for c in module_name) % 1000
    return np.random.RandomState(seed + offset)
```

---

## Reproducibility Checklist

### For Any Analysis

- [ ] Set explicit seed in config (not null)
- [ ] Record Python version
- [ ] Record package versions (`pip freeze`)
- [ ] Archive input data or record hashes
- [ ] Document any preprocessing steps

### For Publication

- [ ] Include seed in methods section
- [ ] Provide complete config file
- [ ] Archive code version (git commit hash)
- [ ] Archive input data (or provide access)
- [ ] Run reproducibility verification

---

## Verification

### Hash Verification

The framework can generate SHA-256 hashes of outputs for verification:

```python
import hashlib

def hash_file(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()

# Compare outputs from two runs
hash1 = hash_file('run_a/pathway_scores.csv')
hash2 = hash_file('run_b/pathway_scores.csv')
assert hash1 == hash2, "Outputs differ!"
```

### Reproducibility Test

Run the same analysis twice and compare:

```bash
# Run twice with same seed
pathway-subtyping --config configs/test_synthetic.yaml
mv outputs/synthetic_test outputs/run_a

pathway-subtyping --config configs/test_synthetic.yaml
mv outputs/synthetic_test outputs/run_b

# Compare outputs
diff outputs/run_a/pathway_scores.csv outputs/run_b/pathway_scores.csv
diff outputs/run_a/subtype_assignments.csv outputs/run_b/subtype_assignments.csv
```

Expected: No differences.

---

## Known Sources of Non-Determinism

Some operations may produce slight variations even with fixed seeds:

| Source | Cause | Mitigation |
|--------|-------|------------|
| Different Python versions | Implementation changes | Pin Python version |
| Different NumPy versions | Algorithm changes | Pin package versions |
| Parallel execution | Thread ordering | Use single thread |
| Float accumulation | Order of operations | Accept small epsilon |
| Hash randomization | Python 3 default | Set PYTHONHASHSEED |

### Environment Variables for Full Determinism

```bash
export PYTHONHASHSEED=42
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

---

## Version Pinning

### Creating a Lock File

```bash
# After testing, freeze versions
pip freeze > requirements.lock
```

### Using a Lock File

```bash
pip install -r requirements.lock
```

### Key Dependencies to Pin

| Package | Purpose |
|---------|---------|
| numpy | Array operations |
| pandas | Data frames |
| scikit-learn | GMM clustering |
| scipy | Statistical functions |

---

## Cross-Platform Considerations

### Potential Differences

- Windows vs. Linux vs. macOS may have subtle differences in:
  - Floating-point arithmetic
  - File path handling
  - Random number generation

### Mitigation

1. Document the platform used
2. Use Docker for exact environment replication
3. Accept small numerical differences (< 1e-10)

### Docker for Reproducibility

```dockerfile
FROM python:3.11-slim

WORKDIR /app
COPY requirements.lock .
RUN pip install -r requirements.lock

COPY . .
RUN pip install -e .

ENV PYTHONHASHSEED=42
CMD ["pathway-subtyping", "--config", "configs/test_synthetic.yaml"]
```

---

## Recording Metadata

### Automatic Metadata

The framework records run metadata in `report.json`:

```json
{
  "run_name": "my_analysis",
  "timestamp": "2026-01-29T10:30:00",
  "seed": 42,
  "summary": {
    "n_samples": 60,
    "n_pathways": 4,
    "n_clusters": 4
  }
}
```

### Manual Metadata

For full reproducibility, also record:

```yaml
# run_metadata.yaml
environment:
  python_version: "3.11.0"
  os: "macOS 14.0"

packages:
  numpy: "1.26.0"
  pandas: "2.1.0"
  scikit-learn: "1.3.0"

git:
  commit: "abc123def456"
  branch: "main"
  clean: true

input_hashes:
  vcf: "sha256:abc123..."
  pathways: "sha256:def456..."
```

---

## Troubleshooting

### Results Differ Between Runs

1. **Check seed is set**:
   ```yaml
   pipeline:
     seed: 42  # Must be integer, not null
   ```

2. **Clear cached data**:
   ```bash
   rm -rf outputs/<run_name>
   ```

3. **Check for external dependencies**:
   - Database connections may change
   - Network resources may differ

### Results Differ Between Machines

1. **Check Python version**:
   ```bash
   python --version
   ```

2. **Check package versions**:
   ```bash
   pip freeze | grep -E "numpy|pandas|scikit-learn"
   ```

3. **Set environment variables**:
   ```bash
   export PYTHONHASHSEED=42
   export OMP_NUM_THREADS=1
   ```

### Small Numerical Differences

Small differences (< 1e-6) in floating-point values are normal and acceptable due to:
- Order of operations in parallel computations
- Platform-specific optimizations
- Compiler differences

Use approximate comparison:
```python
import numpy as np
np.allclose(result1, result2, rtol=1e-5, atol=1e-8)
```

---

## Best Practices Summary

1. **Always set a seed** for reproducible research
2. **Pin package versions** using requirements.lock
3. **Record environment** details with each run
4. **Verify reproducibility** by running twice
5. **Archive inputs** or record their hashes
6. **Document preprocessing** steps
7. **Use version control** for code and configs

---

## See Also

- [Quickstart Guide](quickstart.md) - Getting started
- [Configuration](api/config.md) - Config file reference
- [Troubleshooting](troubleshooting.md) - Common issues
