# CLI API Reference

The command-line interface provides easy access to the pipeline.

## Commands

### `psf` / `pathway-subtyping`

Both commands are aliases for the same CLI.

```bash
psf --config configs/my_config.yaml
pathway-subtyping --config configs/my_config.yaml
```

---

## Options

### `--config PATH`

Path to YAML configuration file (required).

```bash
psf --config configs/test_synthetic.yaml
```

### `--version`

Display version information.

```bash
psf --version
# pathway-subtyping 0.1.0
```

### `--help`

Display help message.

```bash
psf --help
```

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (see stderr for details) |

---

## Examples

### Basic Usage

```bash
# Run with sample data
psf --config configs/test_synthetic.yaml

# View outputs
cat outputs/synthetic_test/report.md
```

### Custom Configuration

```bash
# Copy and edit config
cp configs/example_autism.yaml configs/my_analysis.yaml
# Edit my_analysis.yaml with your paths

# Run
psf --config configs/my_analysis.yaml
```

### Batch Processing

```bash
# Run multiple cohorts
for config in configs/cohort_*.yaml; do
    echo "Running $config..."
    psf --config "$config"
done
```

### With Output Redirection

```bash
# Save logs
psf --config my_config.yaml > pipeline.log 2>&1

# Check for errors
grep -i error pipeline.log
```

---

## Environment Variables

| Variable | Description |
|----------|-------------|
| `PYTHONPATH` | Add `src/` if running from source |

### Running from Source

```bash
# From repository root
export PYTHONPATH=src
python -m pathway_subtyping --config configs/test_synthetic.yaml
```

---

## Module Execution

The CLI can also be invoked as a Python module:

```bash
python -m pathway_subtyping --config configs/test_synthetic.yaml
```
