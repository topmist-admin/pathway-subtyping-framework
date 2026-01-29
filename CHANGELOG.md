# Changelog

All notable changes to the Pathway Subtyping Framework will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- PyPI package publishing preparation

## [0.1.0] - 2026-01-29

### Added

#### Core Pipeline
- Complete pathway subtyping pipeline with VCF → clustering → report workflow
- Gene burden computation with LoF/missense weighting and CADD score normalization
- Pathway score aggregation using GMT file definitions
- GMM clustering with automatic cluster selection via BIC
- Cluster labeling based on top contributing pathways

#### Validation Gates
- Negative Control 1: Label shuffle test (ARI should be < 0.15)
- Negative Control 2: Random gene sets test (ARI should be < 0.15)
- Stability Test: Bootstrap resampling (ARI should be >= 0.8)
- Comprehensive validation reporting in JSON and Markdown

#### Disease Support
- Autism pathway definitions (15 pathways, ~200 genes) - Validated
- Schizophrenia pathway template
- Epilepsy pathway template
- Intellectual disability pathway template
- Guide for adapting to new diseases

#### Infrastructure
- `pyproject.toml` for modern Python packaging
- CLI entry points (`psf`, `pathway-subtyping`)
- YAML-based configuration system
- Reproducibility features (seed control, metadata logging)

#### Testing
- 64 unit and integration tests
- Test fixtures for synthetic data
- CI/CD with GitHub Actions
- Multi-OS (Ubuntu, macOS) and multi-Python (3.9-3.12) testing

#### Documentation
- Comprehensive README with quick start guide
- Disease adaptation guide
- Pathway curation guide
- Validation gates explanation
- Getting-started Jupyter notebook tutorial

#### Containerization
- Multi-stage Dockerfile (runtime, dev, jupyter targets)
- docker-compose.yml for easy orchestration
- Pre-configured development environment

#### Sample Data
- Synthetic VCF with 60 samples and 30 variants
- Synthetic phenotypes with 4 planted subtypes
- Test configuration ready to run out of the box

### Technical Details
- Python 3.8+ support
- Dependencies: numpy, pandas, scikit-learn, scipy, pysam, matplotlib, seaborn
- MIT License

## [0.0.1] - 2026-01-29

### Added
- Initial project structure
- Core module scaffolding
- Basic documentation

---

## Version History Summary

| Version | Date | Highlights |
|---------|------|------------|
| 0.1.0 | 2026-01-29 | First public release with full pipeline |
| 0.0.1 | 2026-01-29 | Initial project setup |

[Unreleased]: https://github.com/topmist-admin/pathway-subtyping-framework/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/topmist-admin/pathway-subtyping-framework/releases/tag/v0.1.0
[0.0.1]: https://github.com/topmist-admin/pathway-subtyping-framework/releases/tag/v0.0.1
