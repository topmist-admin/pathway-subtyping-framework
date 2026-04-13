# Pathway Subtyping Framework
# Multi-stage Docker build for reproducible analysis environments
#
# Stages:
#   runtime     — Minimal production image (CLI + core + viz + graph)
#   development — Adds tests, linters, notebooks
#   jupyter     — Notebook server on port 8888
#
# Build examples:
#   docker build --target runtime -t psf:runtime .
#   docker build --target jupyter -t psf:jupyter .

# =============================================================================
# Stage 1: Builder — compile dependencies with build tools
# =============================================================================
FROM python:3.11-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy only what pip needs to resolve and install
COPY pyproject.toml setup.py README.md LICENSE CITATION.cff ./
COPY src/ ./src/

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir ".[viz,graph]"

# =============================================================================
# Stage 2: Runtime — minimal production image
# =============================================================================
FROM python:3.11-slim AS runtime

LABEL maintainer="Rohit Chauhan <info@topmist.com>"
LABEL description="Pathway Subtyping Framework — disease-agnostic molecular subtype discovery"
LABEL version="0.5.0"
LABEL org.opencontainers.image.source="https://codeberg.org/pathways/pathway-subtyping-framework"
LABEL org.opencontainers.image.documentation="RRID:SCR_028051"

# Runtime C-library dependencies (no -dev packages)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-1.0 \
    liblzma5 \
    && rm -rf /var/lib/apt/lists/*

# Copy pre-built virtualenv from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Create non-root user
RUN useradd --create-home --shell /bin/bash psf
USER psf
WORKDIR /home/psf

# Copy configs (package already installed in venv from builder stage)
COPY --chown=psf:psf configs/ ./configs/

RUN mkdir -p outputs

ENTRYPOINT ["python", "-m", "pathway_subtyping"]
CMD ["--help"]

# =============================================================================
# Stage 3: Development — tests, linters, examples
# =============================================================================
FROM runtime AS development

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Pin linter versions to match CI
RUN pip install --no-cache-dir \
    pytest \
    pytest-cov \
    "black>=26.1,<27" \
    "isort>=7.0,<8" \
    "flake8>=7.3,<8" \
    jupyter \
    ipykernel \
    nbconvert

USER psf

COPY --chown=psf:psf tests/ ./tests/
COPY --chown=psf:psf data/ ./data/
COPY --chown=psf:psf examples/ ./examples/

ENTRYPOINT []
CMD ["/bin/bash"]

# =============================================================================
# Stage 4: Jupyter — notebook server for interactive analysis
# =============================================================================
FROM development AS jupyter

USER root
RUN pip install --no-cache-dir lifelines
USER psf

EXPOSE 8888

ENTRYPOINT []
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--NotebookApp.token=''"]
