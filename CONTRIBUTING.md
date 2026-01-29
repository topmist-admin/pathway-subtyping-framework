# Contributing to Pathway Subtyping Framework

Thank you for your interest in contributing! This document provides guidelines for contributing to the project.

## Ways to Contribute

| Contribution Type | Description |
|-------------------|-------------|
| **Bug Reports** | Report issues you encounter |
| **Feature Requests** | Suggest new features |
| **Pathway Definitions** | Add disease-specific pathway files |
| **Code Contributions** | Fix bugs or add features |
| **Documentation** | Improve docs, tutorials, examples |
| **Testing** | Add tests, report test failures |

## Getting Started

### 1. Fork and Clone

```bash
# Fork via GitHub UI, then:
git clone https://github.com/YOUR_USERNAME/pathway-subtyping-framework
cd pathway-subtyping-framework
git remote add upstream https://github.com/topmist-admin/pathway-subtyping-framework
```

### 2. Set Up Development Environment

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # or venv\Scripts\activate on Windows

# Install with dev dependencies
pip install -e ".[dev]"

# Set up pre-commit hooks
pre-commit install
```

### 3. Create a Branch

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/your-bug-fix
```

## Development Workflow

### Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_validation.py -v

# Run with coverage
pytest tests/ --cov=src/pathway_subtyping --cov-report=html

# Quick smoke test
pytest tests/test_integration.py::TestEndToEndPipeline::test_synthetic_data_pipeline -v
```

### Code Style

We use automated formatting and linting:

```bash
# Format code
black src/ tests/
isort src/ tests/

# Check linting
flake8 src/ tests/

# Type checking
mypy src/
```

Pre-commit hooks will run these automatically on commit.

### Running the Pipeline

```bash
# Run with test data
PYTHONPATH=src python -m pathway_subtyping --config configs/test_synthetic.yaml

# Or after installation
psf --config configs/test_synthetic.yaml
```

## Contribution Guidelines

### Code Contributions

1. **Write tests** for new functionality
2. **Follow existing patterns** in the codebase
3. **Document public APIs** with docstrings
4. **Keep changes focused** — one feature/fix per PR
5. **Update CHANGELOG.md** for user-facing changes

### Code Style Requirements

- Follow PEP 8
- Use type hints for function signatures
- Write docstrings for public functions (Google style)
- Maximum line length: 100 characters
- Keep functions focused and small

### Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
feat: add cross-cohort validation utilities
fix: handle empty pathway scores gracefully
docs: add API documentation for ValidationGates
test: add integration tests for large cohorts
refactor: simplify gene burden computation
chore: update dependencies
```

### Pull Request Process

1. **Update your branch** with latest upstream:
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. **Ensure tests pass**:
   ```bash
   pytest tests/ -v
   ```

3. **Push your branch**:
   ```bash
   git push origin feature/your-feature-name
   ```

4. **Open PR** via GitHub with:
   - Clear title following conventional commits
   - Description of changes
   - Link to related issues
   - Screenshots/outputs for UI changes

5. **Address review feedback** promptly

6. **Squash and merge** when approved

## Contributing Pathway Definitions

See [Contributing Pathways Guide](docs/guides/contributing-pathways.md) for detailed instructions.

### Quick Checklist

- [ ] GMT format with header comments
- [ ] 8+ pathways with 10+ genes each
- [ ] Sources documented (publications, databases)
- [ ] Gene symbols are HGNC standard
- [ ] Validated with GMT parser
- [ ] Config file created (`configs/example_[disease].yaml`)
- [ ] README.md table updated

### Adding a New Disease

1. Create `data/pathways/[disease]_pathways.gmt`
2. Create `configs/example_[disease].yaml`
3. Add entry to README.md table
4. (Optional) Add disease-specific notes in docs

## Reporting Issues

### Bug Reports

Include:
- **Environment**: OS, Python version, package version
- **Steps to reproduce**: Minimal example
- **Expected behavior**: What should happen
- **Actual behavior**: What happens instead
- **Error messages**: Full traceback
- **Data samples**: If possible, synthetic data reproducing the issue

### Feature Requests

Include:
- **Use case**: Why is this needed?
- **Proposed solution**: How should it work?
- **Alternatives considered**: Other approaches

## Testing Guidelines

### Test Categories

| Category | Location | Purpose |
|----------|----------|---------|
| Unit tests | `tests/test_*.py` | Test individual functions |
| Integration | `tests/test_integration.py` | End-to-end pipeline tests |
| Config tests | `tests/test_config.py` | Configuration loading/validation |

### Writing Tests

```python
import pytest
from pathway_subtyping import DemoPipeline, PipelineConfig

def test_my_feature():
    """Test description."""
    # Arrange
    config = PipelineConfig(...)

    # Act
    result = some_function()

    # Assert
    assert result == expected
```

### Test Fixtures

Common fixtures are in `tests/conftest.py`:
- `sample_pathway_scores`: Sample pathway score DataFrame
- `sample_cluster_labels`: Sample cluster label array
- `sample_pathways`: Sample pathway dictionary

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Focus on constructive feedback
- Assume good intentions
- Keep discussions technical and productive

## Questions?

- **GitHub Issues**: For bugs and feature requests
- **GitHub Discussions**: For questions and ideas
- **Email**: info@topmist.com

## Recognition

Contributors are recognized in:
- CONTRIBUTORS.md file
- Release notes for significant contributions
- Publication acknowledgments for major contributions

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for helping improve pathway subtyping research!
