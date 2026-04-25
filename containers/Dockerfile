# =============================================================================
#  Genetic QC Pipeline — Docker image
# =============================================================================
#  Contains every tool required by both the SNP array and WGS/WES QC pipelines.
#  Build:
#    docker build -t genetic-qc:1.0 -f containers/Dockerfile .
#
#  Test:
#    docker run --rm genetic-qc:1.0 plink   --version
#    docker run --rm genetic-qc:1.0 plink2  --version
#    docker run --rm genetic-qc:1.0 samtools --version
#    docker run --rm genetic-qc:1.0 bcftools --version
#    docker run --rm genetic-qc:1.0 gatk     --version
#    docker run --rm genetic-qc:1.0 picard   --version
#    docker run --rm genetic-qc:1.0 fastqc   --version
#    docker run --rm genetic-qc:1.0 mosdepth --version
#    docker run --rm genetic-qc:1.0 python3  --version
#    docker run --rm genetic-qc:1.0 Rscript  --version
#
#  Push to registry (optional):
#    docker tag genetic-qc:1.0 ghcr.io/kaiyao28/genetic-qc:1.0
#    docker push ghcr.io/kaiyao28/genetic-qc:1.0
#
#  Convert to Singularity SIF for HPC:
#    apptainer build genetic-qc.sif docker://ghcr.io/kaiyao28/genetic-qc:1.0
# =============================================================================

FROM condaforge/mambaforge:23.11.0-0

LABEL maintainer="genetic_qc_pipeline" \
      description="All tools for SNP array and WGS/WES genetic QC" \
      version="1.0"

# Prevent interactive prompts during package install
ARG DEBIAN_FRONTEND=noninteractive

# ── System dependencies ───────────────────────────────────────────────────────
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        procps \
        curl \
        wget \
        bzip2 \
        ca-certificates \
        libgomp1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── Core bioinformatics tools via conda/mamba ─────────────────────────────────
# Versions are pinned for reproducibility.
# Update version strings here when upgrading the image.
RUN mamba install -y \
        -c conda-forge \
        -c bioconda \
        "plink=1.90b6.21" \
        "plink2=2.00a5.12" \
        "samtools=1.18" \
        "bcftools=1.18" \
        "picard=3.1.1" \
        "fastqc=0.12.1" \
        "mosdepth=0.3.6" \
        "gatk4=4.5.0.0" \
        "tabix=1.11" \
        "htslib=1.18" \
    && mamba clean -a -y

# ── VerifyBamID2 (contamination check) ───────────────────────────────────────
# Installed separately because it occasionally has dependency conflicts.
# Pipeline falls back to GATK CalculateContamination if not available.
RUN mamba install -y \
        -c conda-forge \
        -c bioconda \
        "verifybamid2=2.0.1" \
    && mamba clean -a -y \
    || echo "WARNING: VerifyBamID2 install failed — contamination step will use GATK fallback"

# ── Python environment ────────────────────────────────────────────────────────
RUN mamba install -y \
        -c conda-forge \
        "python=3.10" \
        "pandas=2.0" \
        "numpy=1.26" \
        "scipy=1.11" \
        "matplotlib=3.8" \
        "seaborn=0.13" \
    && mamba clean -a -y

# ── R environment ─────────────────────────────────────────────────────────────
RUN mamba install -y \
        -c conda-forge \
        "r-base=4.3" \
        "r-ggplot2=3.4" \
        "r-data.table=1.14" \
        "r-dplyr=1.1" \
        "r-readr=2.1" \
        "r-scales=1.3" \
    && mamba clean -a -y

# ── Verify key tool installations ─────────────────────────────────────────────
RUN plink   --version 2>&1 | head -1 && \
    plink2  --version 2>&1 | head -1 && \
    samtools --version 2>&1 | head -1 && \
    bcftools --version 2>&1 | head -1 && \
    gatk --version 2>&1 | head -1 && \
    fastqc  --version 2>&1 | head -1 && \
    mosdepth --version 2>&1 | head -1 && \
    python3 --version 2>&1 && \
    Rscript --version 2>&1

WORKDIR /pipeline

CMD ["/bin/bash"]
