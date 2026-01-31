# Zenodo DOI Setup Guide

This guide walks through the process of minting a DOI for the Pathway Subtyping Framework via Zenodo integration with GitHub.

## Overview

Zenodo is a general-purpose open repository that:
- Provides **persistent DOIs** for research software
- Automatically archives **GitHub releases**
- Is recognized by journals and funders
- Is **free** for open-source projects

## Prerequisites

- GitHub repository (public)
- Zenodo account (free, can use GitHub OAuth)
- At least one GitHub release

---

## Step 1: Create Zenodo Account

1. Go to [https://zenodo.org](https://zenodo.org)
2. Click **Sign Up** (top right)
3. Choose **Sign up with GitHub** for easy integration
4. Authorize Zenodo to access your GitHub account

---

## Step 2: Link GitHub Repository

1. Go to [https://zenodo.org/account/settings/github/](https://zenodo.org/account/settings/github/)
2. Find `topmist-admin/pathway-subtyping-framework` in the repository list
3. Toggle the switch to **ON**

If the repository doesn't appear:
- Click **Sync now** to refresh the list
- Ensure the repository is public
- Check that you have admin access to the repo

---

## Step 3: Configure Zenodo Metadata

Create a `.zenodo.json` file in the repository root:

```json
{
  "title": "Pathway Subtyping Framework",
  "description": "A disease-agnostic framework for pathway-based molecular subtype discovery in genetically heterogeneous conditions.",
  "creators": [
    {
      "name": "Chauhan, Rohit",
      "affiliation": "Independent Researcher",
      "orcid": "0000-0000-0000-0000"
    }
  ],
  "keywords": [
    "bioinformatics",
    "genomics",
    "pathway analysis",
    "molecular subtyping",
    "rare variants",
    "clustering",
    "precision medicine",
    "genetic heterogeneity"
  ],
  "license": "MIT",
  "upload_type": "software",
  "access_right": "open",
  "related_identifiers": [
    {
      "identifier": "https://github.com/topmist-admin/pathway-subtyping-framework",
      "relation": "isSupplementTo",
      "scheme": "url"
    }
  ],
  "communities": [
    {"identifier": "bioinformatics"},
    {"identifier": "genomics"}
  ]
}
```

### Metadata Fields

| Field | Description | Required |
|-------|-------------|----------|
| title | Project name | Yes |
| description | Brief summary | Yes |
| creators | Author list with affiliations | Yes |
| keywords | Search terms | Recommended |
| license | License identifier | Recommended |
| upload_type | "software" for code | Yes |
| access_right | "open" for public | Yes |
| communities | Zenodo community IDs | Optional |

### Finding Your ORCID

If you don't have an ORCID:
1. Go to [https://orcid.org/register](https://orcid.org/register)
2. Create a free account
3. Add your ORCID to the `.zenodo.json`

---

## Step 4: Create GitHub Release

### Via GitHub Web Interface

1. Go to repository → **Releases** → **Create a new release**
2. Click **Choose a tag** → type `v0.1.0` → **Create new tag**
3. Set **Release title**: `v0.1.0 - Initial Release`
4. Add **Release notes**:

```markdown
## Pathway Subtyping Framework v0.1.0

Initial release of the disease-agnostic pathway subtyping framework.

### Features
- VCF parsing with gene/consequence annotation
- Pathway-level burden scoring
- GMM-based clustering with BIC selection
- Three validation gates (label shuffle, random genes, bootstrap)
- Cross-cohort comparison utilities
- CLI and Python API

### Installation
```bash
pip install pathway-subtyping
```

### Documentation
- [Quickstart Guide](docs/quickstart.md)
- [API Reference](docs/api/index.md)

### Requirements
- Python 3.9+
- See requirements.txt for dependencies
```

5. Click **Publish release**

### Via Command Line

```bash
# Create and push tag
git tag -a v0.1.0 -m "Initial release"
git push origin v0.1.0

# Create release via GitHub CLI
gh release create v0.1.0 \
  --title "v0.1.0 - Initial Release" \
  --notes-file RELEASE_NOTES.md
```

---

## Step 5: Verify DOI Minting

After creating the GitHub release:

1. Wait 5-10 minutes for Zenodo to process
2. Go to [https://zenodo.org/account/settings/github/](https://zenodo.org/account/settings/github/)
3. Find your repository in the list
4. Click the DOI badge to view the record

The DOI will look like: `10.5281/zenodo.XXXXXXX`

---

## Step 6: Add DOI Badge to README

Once you have the DOI, add badges to your README.md:

```markdown
# Pathway Subtyping Framework

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![PyPI version](https://badge.fury.io/py/pathway-subtyping.svg)](https://badge.fury.io/py/pathway-subtyping)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
```

Replace `XXXXXXX` with your actual Zenodo record number.

---

## Step 7: Citation File

Create a `CITATION.cff` file for standardized citation:

```yaml
cff-version: 1.2.0
title: Pathway Subtyping Framework
message: >-
  If you use this software, please cite it using
  the metadata from this file.
type: software
authors:
  - family-names: Chauhan
    given-names: Rohit
    orcid: 'https://orcid.org/0000-0000-0000-0000'
repository-code: 'https://github.com/topmist-admin/pathway-subtyping-framework'
url: 'https://github.com/topmist-admin/pathway-subtyping-framework'
license: MIT
version: 0.1.0
date-released: '2026-01-29'
doi: '10.5281/zenodo.XXXXXXX'
keywords:
  - bioinformatics
  - genomics
  - pathway analysis
  - molecular subtyping
```

---

## Versioning Strategy

### Semantic Versioning

Follow [SemVer](https://semver.org/):
- **MAJOR** (1.0.0): Breaking API changes
- **MINOR** (0.2.0): New features, backward compatible
- **PATCH** (0.1.1): Bug fixes, backward compatible

### DOI Strategy

Zenodo provides two types of DOIs:

1. **Version-specific DOI**: Points to exact version (e.g., `10.5281/zenodo.1234567`)
2. **Concept DOI**: Points to all versions (e.g., `10.5281/zenodo.1234566`)

Use the **concept DOI** in citations for longevity.

---

## Troubleshooting

### Release Not Appearing in Zenodo

1. Check webhook: Repository Settings → Webhooks → Zenodo
2. Verify the repository toggle is ON in Zenodo settings
3. Check Zenodo's GitHub integration page for errors
4. Wait 10-15 minutes and refresh

### DOI Not Resolving

- New DOIs may take up to 24 hours to resolve
- Check the Zenodo record page directly first

### Updating Metadata After Release

1. Go to your Zenodo record
2. Click **Edit**
3. Update metadata
4. Click **Save** and **Publish**

Note: This creates a new version of the record but keeps the same concept DOI.

---

## Best Practices

1. **Release early, release often**: Each release gets a DOI
2. **Write good release notes**: They appear in Zenodo
3. **Use concept DOI in papers**: Points to latest version
4. **Keep `.zenodo.json` updated**: Metadata propagates to new releases
5. **Add ORCID**: Links your work to your identity

---

## Resources

- [Zenodo Documentation](https://help.zenodo.org/)
- [GitHub-Zenodo Integration](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content)
- [Semantic Versioning](https://semver.org/)
- [CITATION.cff Format](https://citation-file-format.github.io/)

---

## Quick Reference

| Task | Command/Action |
|------|----------------|
| Create tag | `git tag -a v0.1.0 -m "Release v0.1.0"` |
| Push tag | `git push origin v0.1.0` |
| Create release | `gh release create v0.1.0` |
| Check DOI | Visit Zenodo GitHub settings |

---

*Last updated: January 2026*
