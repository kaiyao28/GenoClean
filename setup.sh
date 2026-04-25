#!/usr/bin/env bash
# =============================================================================
#  setup.sh — install all pipeline dependencies
# =============================================================================
#  Usage:
#    bash setup.sh              # install via Mamba (default, works everywhere)
#    bash setup.sh docker       # build Docker image instead
#    bash setup.sh singularity  # build Docker image then convert to SIF
#
#  After setup, run:  bash test_env.sh
# =============================================================================
set -euo pipefail

MODE="${1:-conda}"
IMAGE="genetic-qc:1.0"
SIF="containers/genetic-qc.sif"
ENV_FILE="containers/environment.yml"
DOCKERFILE="containers/Dockerfile"

# ── helpers ───────────────────────────────────────────────────────────────────
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[OK]${NC}  $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
err()  { echo -e "${RED}[ERR]${NC} $*"; exit 1; }

echo ""
echo "=== Genetic QC Pipeline — Setup (mode: ${MODE}) ==="
echo ""

# ── Check Nextflow ─────────────────────────────────────────────────────────────
if ! command -v nextflow &>/dev/null; then
    warn "Nextflow not found — installing..."
    curl -s https://get.nextflow.io | bash
    mkdir -p "$HOME/bin"
    mv nextflow "$HOME/bin/"
    export PATH="$HOME/bin:$PATH"
    ok "Nextflow installed: $(nextflow -version 2>&1 | head -1)"
else
    ok "Nextflow: $(nextflow -version 2>&1 | head -1)"
fi

# ── Mode: conda (default) ──────────────────────────────────────────────────────
if [[ "$MODE" == "conda" ]]; then
    # Prefer mamba, fall back to conda
    PKG_MGR=""
    command -v mamba &>/dev/null && PKG_MGR="mamba"
    command -v conda &>/dev/null && PKG_MGR="${PKG_MGR:-conda}"
    [[ -z "$PKG_MGR" ]] && err "Neither mamba nor conda found. Install Miniforge: https://github.com/conda-forge/miniforge"

    ok "Using: $PKG_MGR"

    if $PKG_MGR env list | grep -q "^genetic_qc "; then
        warn "Environment 'genetic_qc' already exists — updating..."
        $PKG_MGR env update -n genetic_qc -f "$ENV_FILE" --prune
    else
        echo "Creating environment from $ENV_FILE ..."
        $PKG_MGR env create -f "$ENV_FILE"
    fi

    ok "Environment ready. Activate with:  mamba activate genetic_qc"
    echo ""
    echo "Then run the pipeline with:  -profile conda"

# ── Mode: docker ──────────────────────────────────────────────────────────────
elif [[ "$MODE" == "docker" ]]; then
    command -v docker &>/dev/null || err "Docker not found. Install from https://docs.docker.com/engine/install/"
    echo "Building Docker image ${IMAGE} ..."
    docker build -t "$IMAGE" -f "$DOCKERFILE" .
    ok "Docker image built: ${IMAGE}"
    echo ""
    echo "Run the pipeline with:  -profile docker"

# ── Mode: singularity ─────────────────────────────────────────────────────────
elif [[ "$MODE" == "singularity" ]]; then
    # Build Docker first, then convert
    APPTAINER=""
    command -v apptainer   &>/dev/null && APPTAINER="apptainer"
    command -v singularity &>/dev/null && APPTAINER="${APPTAINER:-singularity}"
    [[ -z "$APPTAINER" ]] && err "Apptainer/Singularity not found."

    command -v docker &>/dev/null || err "Docker is required to build the image before converting to SIF."

    echo "Building Docker image ${IMAGE} ..."
    docker build -t "$IMAGE" -f "$DOCKERFILE" .
    ok "Docker image built"

    echo "Converting to SIF: ${SIF} ..."
    mkdir -p containers
    $APPTAINER build "$SIF" "docker-daemon://${IMAGE}"
    ok "SIF built: ${SIF}"
    echo ""
    echo "Run the pipeline with:  -profile singularity"

else
    err "Unknown mode '${MODE}'. Use: conda | docker | singularity"
fi

echo ""
echo "=== Setup complete. Now run:  bash test_env.sh ==="
echo ""
