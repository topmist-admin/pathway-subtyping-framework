# Contributing to Pathway Subtyping Framework

Thank you for your interest in contributing! This document provides guidelines for contributing to the project.

## Ways to Contribute

### 1. Add Disease Pathway Definitions
We welcome pathway GMT files for additional diseases:
- Curate genes from literature and databases
- Follow the [Pathway Curation Guide](docs/guides/pathway-curation-guide.md)
- Submit a PR with your GMT file and documentation

### 2. Improve Documentation
- Fix typos or unclear explanations
- Add examples for your disease area
- Translate documentation

### 3. Report Issues
- Bug reports with reproducible examples
- Feature requests with use case descriptions
- Questions about usage

### 4. Code Contributions
- Performance optimizations
- New features (discuss first via issue)
- Bug fixes

## Development Setup

```bash
# Clone the repository
git clone https://github.com/topmist-admin/pathway-subtyping-framework
cd pathway-subtyping-framework

# Create virtual environment
python -m venv venv
source venv/bin/activate  # or `venv\Scripts\activate` on Windows

# Install in development mode
pip install -r requirements.txt
pip install -e ".[dev]"

# Run tests
pytest tests/
```

## Code Style

- Follow PEP 8
- Use type hints
- Write docstrings for public functions
- Keep functions focused and small

## Pull Request Process

1. **Fork** the repository
2. **Create a branch** for your feature (`git checkout -b feat/my-feature`)
3. **Make changes** with clear commits
4. **Write tests** for new functionality
5. **Update documentation** if needed
6. **Submit PR** with description of changes

### Commit Messages

Use conventional commits:
```
feat: add parkinson's pathway definitions
fix: handle missing CADD scores gracefully
docs: clarify GMT file format
test: add bootstrap stability edge cases
```

## Adding a New Disease

To add pathway definitions for a new disease:

1. Create `data/pathways/[disease]_pathways.gmt`
2. Create `configs/example_[disease].yaml`
3. Add entry to README.md table
4. (Optional) Add to docs with disease-specific notes

### GMT File Requirements
- Tab-separated format
- Gene symbols must be HGNC standard
- Include source URL in description field
- 5-50 genes per pathway recommended

## Testing

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=pathway_subtyping

# Run specific test file
pytest tests/test_clustering.py
```

## Questions?

- Open a [GitHub Discussion](https://github.com/topmist-admin/pathway-subtyping-framework/discussions)
- Email: info@topmist.com

## Code of Conduct

Be respectful and constructive. We're all here to advance research.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
