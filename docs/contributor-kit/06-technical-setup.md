# Technical Setup Guide

This guide covers setting up your development environment for the Pathway Subtyping Framework.

## System Requirements

### Minimum Requirements
- **OS**: macOS 10.15+, Ubuntu 20.04+, Windows 10+ (WSL2)
- **Python**: 3.9 or higher
- **RAM**: 8 GB minimum (16 GB recommended)
- **Storage**: 10 GB free (more for large datasets)

### Recommended for Large Datasets
- **RAM**: 32 GB+
- **CPU**: 8+ cores
- **Storage**: SSD with 100+ GB free

## Local Development Setup

### Step 1: Install Python

**macOS (using Homebrew):**
```bash
brew install python@3.11
```

**Ubuntu:**
```bash
sudo apt update
sudo apt install python3.11 python3.11-venv python3-pip
```

**Windows (WSL2):**
```bash
sudo apt update
sudo apt install python3.11 python3.11-venv python3-pip
```

### Step 2: Clone Repository

```bash
git clone https://github.com/topmist-admin/pathway-subtyping-framework.git
cd pathway-subtyping-framework
```

### Step 3: Create Virtual Environment

```bash
# Create virtual environment
python3.11 -m venv venv

# Activate (macOS/Linux)
source venv/bin/activate

# Activate (Windows)
.\venv\Scripts\activate
```

### Step 4: Install Dependencies

```bash
# Install framework in development mode
pip install -e .

# Or install from requirements
pip install -r requirements.txt
```

### Step 5: Verify Installation

```bash
# Check version
python -c "import pathway_subtyping; print(pathway_subtyping.__version__)"

# Run on synthetic data
psf --config configs/test_synthetic.yaml
```

## IDE Setup

### VS Code (Recommended)

**Install extensions:**
- Python (Microsoft)
- Pylance
- Jupyter
- GitLens

**Workspace settings** (`.vscode/settings.json`):
```json
{
    "python.defaultInterpreterPath": "${workspaceFolder}/venv/bin/python",
    "python.linting.enabled": true,
    "python.linting.pylintEnabled": true,
    "python.formatting.provider": "black",
    "editor.formatOnSave": true,
    "python.testing.pytestEnabled": true
}
```

### PyCharm

1. Open project folder
2. Configure interpreter: Settings → Project → Python Interpreter → Add → Existing environment → Select `venv/bin/python`
3. Enable pytest: Settings → Tools → Python Integrated Tools → Default test runner: pytest

## Cloud Setup

### Google Colab (Quickest Start)

No local setup needed. Open the notebook directly:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/topmist-admin/pathway-subtyping-framework/blob/main/examples/notebooks/01_getting_started.ipynb)

### AWS Setup

**For controlled-access data analysis:**

1. **Launch EC2 instance:**
   - AMI: Ubuntu 22.04 LTS
   - Instance type: r5.xlarge (4 vCPU, 32 GB RAM) or larger
   - Storage: 100 GB EBS (encrypted)

2. **Security group:**
   - SSH (port 22) from your IP only
   - No other inbound traffic

3. **Connect and setup:**
```bash
ssh -i your-key.pem ubuntu@ec2-xxx.compute.amazonaws.com

# Install dependencies
sudo apt update && sudo apt install -y python3.11 python3.11-venv git

# Clone and setup
git clone https://github.com/topmist-admin/pathway-subtyping-framework.git
cd pathway-subtyping-framework
python3.11 -m venv venv
source venv/bin/activate
pip install -e .
```

4. **Data storage:**
   - Use S3 with encryption for data staging
   - Copy to encrypted EBS for analysis
   - Delete after project completion

### Google Cloud Setup

1. **Create VM:**
   - Machine type: n2-highmem-4 (4 vCPU, 32 GB)
   - Boot disk: Ubuntu 22.04, 100 GB SSD
   - Enable encryption

2. **Setup:**
```bash
# SSH to instance
gcloud compute ssh your-instance-name

# Same setup as AWS
sudo apt update && sudo apt install -y python3.11 python3.11-venv git
# ... continue with clone and setup
```

## Data Management

### Directory Structure

```
~/genetic-research/
├── pathway-subtyping-framework/   # Framework code
├── data/                          # Data files (encrypted)
│   ├── raw/                       # Original VCFs
│   ├── processed/                 # Annotated/filtered
│   └── pathways/                  # GMT files
├── outputs/                       # Analysis outputs
│   ├── autism_study/
│   └── schizophrenia_study/
└── configs/                       # Project configs
```

### Handling Large VCF Files

**For files > 10 GB:**

```bash
# Use tabix for indexed access
tabix -p vcf large_cohort.vcf.gz

# Query specific regions
tabix large_cohort.vcf.gz chr1:1000000-2000000

# Split by chromosome
for chr in {1..22} X Y; do
    bcftools view -r chr${chr} large_cohort.vcf.gz -Oz -o split/chr${chr}.vcf.gz
done
```

## Dependency Notes

### Core Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | 1.24+ | Numerical computing |
| pandas | 2.0+ | Data manipulation |
| scikit-learn | 1.3+ | Clustering, metrics |
| scipy | 1.11+ | Statistical functions |
| pyyaml | 6.0+ | Config parsing |
| click | 8.0+ | CLI interface |

### Optional Dependencies

| Package | Purpose | Install |
|---------|---------|---------|
| pysam | VCF parsing (faster) | `pip install pysam` |
| matplotlib | Visualizations | `pip install matplotlib` |
| seaborn | Statistical plots | `pip install seaborn` |
| jupyter | Notebooks | `pip install jupyter` |

### Troubleshooting Dependencies

**NumPy compatibility (especially in Colab):**
```bash
pip install "numpy>=1.24.0,<2.0.0"
```

**pysam installation issues:**
```bash
# macOS
brew install htslib
pip install pysam

# Ubuntu
sudo apt install libhts-dev
pip install pysam
```

## Testing Your Setup

### Run Unit Tests

```bash
# All tests
pytest tests/

# Specific test file
pytest tests/test_pipeline.py

# With coverage
pytest --cov=pathway_subtyping tests/
```

### Run Demo Analysis

```bash
# Synthetic data demo
psf --config configs/test_synthetic.yaml --verbose

# Check outputs
ls -la outputs/synthetic_test/
```

### Verify GPU (if using)

```python
import torch
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"Device: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'}")
```

*Note: GPU is not required for the framework but can accelerate some operations.*

## Environment Variables

Set these for convenience:

```bash
# Add to ~/.bashrc or ~/.zshrc
export PSF_HOME=~/genetic-research/pathway-subtyping-framework
export PSF_DATA=~/genetic-research/data
export PSF_OUTPUT=~/genetic-research/outputs

# Activate env automatically
alias psf-activate='cd $PSF_HOME && source venv/bin/activate'
```

## Next Steps

After completing setup:

1. Run the demo notebook: `examples/notebooks/01_getting_started.ipynb`
2. Review analysis workflow: [07-analysis-workflow.md](07-analysis-workflow.md)
3. Start your data access application: [02-data-access-guide.md](02-data-access-guide.md)

---

*Setup issues? Post in the team channel or create a GitHub issue.*
