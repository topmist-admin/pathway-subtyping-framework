# Pathway Subtyping Framework
# Multi-stage Docker build for reproducible analysis environments

# =============================================================================
# Stage 1: Builder - Install dependencies
# =============================================================================
FROM python:3.11-slim as builder

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Create virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# =============================================================================
# Stage 2: Runtime - Minimal production image
# =============================================================================
FROM python:3.11-slim as runtime

LABEL maintainer="Rohit Chauhan <info@topmist.com>"
LABEL description="Pathway Subtyping Framework - Disease-agnostic molecular subtype discovery"
LABEL version="0.1.0"

# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-1.0 \
    liblzma5 \
    && rm -rf /var/lib/apt/lists/*

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Create non-root user
RUN useradd --create-home --shell /bin/bash psf
USER psf
WORKDIR /home/psf

# Copy application code
COPY --chown=psf:psf src/ ./src/
COPY --chown=psf:psf data/pathways/ ./data/pathways/
COPY --chown=psf:psf configs/ ./configs/

# Install package
COPY --chown=psf:psf pyproject.toml setup.py ./
RUN pip install --no-cache-dir -e .

# Create output directory
RUN mkdir -p outputs

# Default command
ENTRYPOINT ["python", "-m", "pathway_subtyping"]
CMD ["--help"]

# =============================================================================
# Stage 3: Development - Full development environment
# =============================================================================
FROM runtime as development

USER root

# Install development dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Install dev Python packages
RUN pip install --no-cache-dir \
    pytest \
    pytest-cov \
    black \
    isort \
    flake8 \
    jupyter \
    ipykernel

USER psf

# Copy test files
COPY --chown=psf:psf tests/ ./tests/
COPY --chown=psf:psf data/sample/ ./data/sample/
COPY --chown=psf:psf examples/ ./examples/

# Default to bash for development
CMD ["/bin/bash"]

# =============================================================================
# Stage 4: Jupyter - Notebook server for interactive analysis
# =============================================================================
FROM development as jupyter

USER psf

# Expose Jupyter port
EXPOSE 8888

# Start Jupyter notebook
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--NotebookApp.token=''"]
