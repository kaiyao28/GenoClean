#!/usr/bin/env bash
# =============================================================================
#  setup_hpc_manual.sh — download pipeline tools without conda/docker
# =============================================================================
#  Use this when conda, Docker, and Singularity are all unavailable (e.g.
#  restricted HPC nodes with no package manager or internet-restricted login).
#
#  What it does:
#    - Downloads each tool as a pre-compiled binary or compiles from source
#    - Creates a single $TOOL_DIR/bin/ directory with all executables
#    - Prints the nextflow command-line flag to point the pipeline at them
#
#  Requirements (usually already on HPC login nodes):
#    - wget or curl
#    - Java 17+ (for Nextflow, GATK, Picard, FastQC)
#    - gcc + make (for samtools, bcftools, htslib compilation)
#      If not available, load them first: module load gcc
#
#  Usage:
#    bash setup_hpc_manual.sh                        # installs to ~/genetic_qc_tools
#    bash setup_hpc_manual.sh --dir /shared/tools    # custom location
#    bash setup_hpc_manual.sh --dir /shared/tools --skip-compile   # skip samtools build
#
#  After running, add this to your nextflow command:
#    --tool_dir $TOOL_DIR   OR   export GENETIC_QC_TOOL_DIR=$TOOL_DIR
#  And use:
#    -profile manual_paths           # local execution
#    -profile slurm,manual_paths     # SLURM + manual tools
# =============================================================================
set -euo pipefail

# ── Parse arguments ────────────────────────────────────────────────────────────
TOOL_DIR="$HOME/genetic_qc_tools"
SKIP_COMPILE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dir)       TOOL_DIR="$2"; shift 2 ;;
        --skip-compile) SKIP_COMPILE=true; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ── Helpers ────────────────────────────────────────────────────────────────────
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[OK]${NC}  $*"; }
warn() { echo -e "${YELLOW}[SKIP]${NC} $*"; }
err()  { echo -e "${RED}[ERR]${NC} $*"; }

# Use wget if available, otherwise curl
fetch() {
    local url="$1" dest="$2"
    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "$dest" "$url"
    elif command -v curl &>/dev/null; then
        curl -L --progress-bar -o "$dest" "$url"
    else
        echo "ERROR: neither wget nor curl found"; exit 1
    fi
}

mkdir -p "$TOOL_DIR/bin" "$TOOL_DIR/src"

echo ""
echo "=== Genetic QC — Manual Tool Setup ==="
echo "Installing to: $TOOL_DIR"
echo ""

# ── 0. Check Java ──────────────────────────────────────────────────────────────
# Java is required for Nextflow, GATK, Picard, and FastQC.
# On HPC: try  module load java  or  module load openjdk  first.
if ! command -v java &>/dev/null; then
    err "Java not found. On HPC try:  module load java"
    err "Or download a portable JDK: https://adoptium.net/temurin/releases/?version=17"
    err "Extract it and add its bin/ to PATH, then re-run this script."
    exit 1
fi
JAVA_VERSION=$(java -version 2>&1 | head -1)
ok "Java: $JAVA_VERSION"

# ── 1. Nextflow ────────────────────────────────────────────────────────────────
if command -v nextflow &>/dev/null; then
    ok "Nextflow already on PATH: $(nextflow -version 2>&1 | head -1)"
else
    echo "Downloading Nextflow..."
    fetch "https://get.nextflow.io" "$TOOL_DIR/bin/nextflow"
    chmod +x "$TOOL_DIR/bin/nextflow"
    ok "Nextflow installed: $TOOL_DIR/bin/nextflow"
fi

# ── 2. PLINK 1.9 ──────────────────────────────────────────────────────────────
# Pre-compiled static Linux x86_64 binary. No compilation needed.
# Source: https://www.cog-genomics.org/plink/
if [[ -f "$TOOL_DIR/bin/plink" ]]; then
    ok "PLINK 1.9 already present"
else
    echo "Downloading PLINK 1.9..."
    mkdir -p "$TOOL_DIR/src/plink1"
    fetch "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_latest.zip" \
          "$TOOL_DIR/src/plink1/plink1.zip"
    unzip -q "$TOOL_DIR/src/plink1/plink1.zip" -d "$TOOL_DIR/src/plink1/"
    cp "$TOOL_DIR/src/plink1/plink" "$TOOL_DIR/bin/plink"
    chmod +x "$TOOL_DIR/bin/plink"
    ok "PLINK 1.9: $TOOL_DIR/bin/plink"
fi

# ── 3. PLINK 2 ────────────────────────────────────────────────────────────────
# Pre-compiled static Linux x86_64 binary.
# Source: https://www.cog-genomics.org/plink/2.0/
# Note: URL below uses the AVX2 build. If your CPU lacks AVX2, use the
# non-AVX2 build (check: grep -c avx2 /proc/cpuinfo)
if [[ -f "$TOOL_DIR/bin/plink2" ]]; then
    ok "PLINK 2 already present"
else
    echo "Downloading PLINK 2..."
    mkdir -p "$TOOL_DIR/src/plink2"
    fetch "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_latest.zip" \
          "$TOOL_DIR/src/plink2/plink2.zip"
    unzip -q "$TOOL_DIR/src/plink2/plink2.zip" -d "$TOOL_DIR/src/plink2/"
    cp "$TOOL_DIR/src/plink2/plink2" "$TOOL_DIR/bin/plink2"
    chmod +x "$TOOL_DIR/bin/plink2"
    ok "PLINK 2: $TOOL_DIR/bin/plink2"
fi

# ── 4. htslib / samtools / bcftools ───────────────────────────────────────────
# Compiled from source. Requires gcc and make.
# If unavailable: module load gcc   (or module load samtools if pre-built)
#
# Versions are pinned to match the Docker image (see containers/Dockerfile).
HTSLIB_VER="1.18"
SAMTOOLS_VER="1.18"
BCFTOOLS_VER="1.18"

if [[ "$SKIP_COMPILE" == "true" ]]; then
    warn "Skipping samtools/bcftools compilation (--skip-compile). Ensure they are on PATH."
else
    if ! command -v gcc &>/dev/null; then
        err "gcc not found — cannot compile samtools/bcftools."
        err "Try:  module load gcc   then re-run this script."
        err "Or use --skip-compile and ensure samtools/bcftools are on PATH another way."
        exit 1
    fi

    # htslib (required by samtools and bcftools)
    if [[ -f "$TOOL_DIR/bin/bgzip" ]]; then
        ok "htslib already built"
    else
        echo "Downloading and building htslib ${HTSLIB_VER}..."
        fetch "https://github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2" \
              "$TOOL_DIR/src/htslib.tar.bz2"
        tar -xjf "$TOOL_DIR/src/htslib.tar.bz2" -C "$TOOL_DIR/src/"
        pushd "$TOOL_DIR/src/htslib-${HTSLIB_VER}" > /dev/null
        ./configure --prefix="$TOOL_DIR" --disable-libcurl 2>&1 | tail -3
        make -j4 2>&1 | tail -3
        make install 2>&1 | tail -3
        popd > /dev/null
        ok "htslib: $TOOL_DIR/bin/bgzip"
    fi

    # samtools
    if [[ -f "$TOOL_DIR/bin/samtools" ]]; then
        ok "samtools already built"
    else
        echo "Downloading and building samtools ${SAMTOOLS_VER}..."
        fetch "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2" \
              "$TOOL_DIR/src/samtools.tar.bz2"
        tar -xjf "$TOOL_DIR/src/samtools.tar.bz2" -C "$TOOL_DIR/src/"
        pushd "$TOOL_DIR/src/samtools-${SAMTOOLS_VER}" > /dev/null
        ./configure --prefix="$TOOL_DIR" --with-htslib="$TOOL_DIR" 2>&1 | tail -3
        make -j4 2>&1 | tail -3
        make install 2>&1 | tail -3
        popd > /dev/null
        ok "samtools: $TOOL_DIR/bin/samtools"
    fi

    # bcftools
    if [[ -f "$TOOL_DIR/bin/bcftools" ]]; then
        ok "bcftools already built"
    else
        echo "Downloading and building bcftools ${BCFTOOLS_VER}..."
        fetch "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VER}/bcftools-${BCFTOOLS_VER}.tar.bz2" \
              "$TOOL_DIR/src/bcftools.tar.bz2"
        tar -xjf "$TOOL_DIR/src/bcftools.tar.bz2" -C "$TOOL_DIR/src/"
        pushd "$TOOL_DIR/src/bcftools-${BCFTOOLS_VER}" > /dev/null
        ./configure --prefix="$TOOL_DIR" --with-htslib="$TOOL_DIR" 2>&1 | tail -3
        make -j4 2>&1 | tail -3
        make install 2>&1 | tail -3
        popd > /dev/null
        ok "bcftools: $TOOL_DIR/bin/bcftools"
    fi
fi

# ── 5. GATK 4 ─────────────────────────────────────────────────────────────────
# The zip includes a  gatk  shell wrapper that calls the JAR via Java.
# Source: https://github.com/broadinstitute/gatk/releases
GATK_VER="4.5.0.0"
if [[ -f "$TOOL_DIR/bin/gatk" ]]; then
    ok "GATK already present"
else
    echo "Downloading GATK ${GATK_VER}..."
    mkdir -p "$TOOL_DIR/src/gatk"
    fetch "https://github.com/broadinstitute/gatk/releases/download/${GATK_VER}/gatk-${GATK_VER}.zip" \
          "$TOOL_DIR/src/gatk/gatk.zip"
    unzip -q "$TOOL_DIR/src/gatk/gatk.zip" -d "$TOOL_DIR/src/gatk/"
    ln -sf "$TOOL_DIR/src/gatk/gatk-${GATK_VER}/gatk" "$TOOL_DIR/bin/gatk"
    chmod +x "$TOOL_DIR/src/gatk/gatk-${GATK_VER}/gatk"
    ok "GATK: $TOOL_DIR/bin/gatk"
fi

# ── 6. Picard ─────────────────────────────────────────────────────────────────
# Downloaded as a single JAR. A wrapper script is created so the pipeline
# can call  picard MarkDuplicates  just like a regular command.
# Source: https://github.com/broadinstitute/picard/releases
PICARD_VER="3.1.1"
if [[ -f "$TOOL_DIR/bin/picard" ]]; then
    ok "Picard already present"
else
    echo "Downloading Picard ${PICARD_VER}..."
    mkdir -p "$TOOL_DIR/picard"
    fetch "https://github.com/broadinstitute/picard/releases/download/${PICARD_VER}/picard.jar" \
          "$TOOL_DIR/picard/picard.jar"
    cat > "$TOOL_DIR/bin/picard" << EOF
#!/usr/bin/env bash
exec java -jar "$TOOL_DIR/picard/picard.jar" "\$@"
EOF
    chmod +x "$TOOL_DIR/bin/picard"
    ok "Picard: $TOOL_DIR/bin/picard (wrapper → $TOOL_DIR/picard/picard.jar)"
fi

# ── 7. FastQC ─────────────────────────────────────────────────────────────────
# The zip includes a  fastqc  shell script that wraps the Java JAR.
# Source: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
FASTQC_VER="0.12.1"
if [[ -f "$TOOL_DIR/bin/fastqc" ]]; then
    ok "FastQC already present"
else
    echo "Downloading FastQC ${FASTQC_VER}..."
    mkdir -p "$TOOL_DIR/src/fastqc"
    fetch "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VER}.zip" \
          "$TOOL_DIR/src/fastqc/fastqc.zip"
    unzip -q "$TOOL_DIR/src/fastqc/fastqc.zip" -d "$TOOL_DIR/src/fastqc/"
    chmod +x "$TOOL_DIR/src/fastqc/FastQC/fastqc"
    ln -sf "$TOOL_DIR/src/fastqc/FastQC/fastqc" "$TOOL_DIR/bin/fastqc"
    ok "FastQC: $TOOL_DIR/bin/fastqc"
fi

# ── 8. mosdepth ───────────────────────────────────────────────────────────────
# Pre-compiled static Linux binary. No dependencies needed.
# Source: https://github.com/brentp/mosdepth/releases
MOSDEPTH_VER="0.3.6"
if [[ -f "$TOOL_DIR/bin/mosdepth" ]]; then
    ok "mosdepth already present"
else
    echo "Downloading mosdepth ${MOSDEPTH_VER}..."
    fetch "https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VER}/mosdepth" \
          "$TOOL_DIR/bin/mosdepth"
    chmod +x "$TOOL_DIR/bin/mosdepth"
    ok "mosdepth: $TOOL_DIR/bin/mosdepth"
fi

# ── 9. VerifyBamID2 ───────────────────────────────────────────────────────────
# Binary + SVD resource files. The pipeline uses VERIFYBAMID2_SVD_PREFIX
# to find the resource files; that env var is set automatically by
# conf/manual_tools.config based on --tool_dir.
#
# Source: https://github.com/Griffan/VerifyBamID/releases
# If download fails, the pipeline falls back to GATK CalculateContamination.
VERIFYBAMID_VER="2.0.1"
if [[ -f "$TOOL_DIR/bin/VerifyBamID" ]]; then
    ok "VerifyBamID2 already present"
else
    echo "Downloading VerifyBamID2 ${VERIFYBAMID_VER}..."
    mkdir -p "$TOOL_DIR/src/verifybamid2" "$TOOL_DIR/verifybamid2"
    VBID_URL="https://github.com/Griffan/VerifyBamID/releases/download/${VERIFYBAMID_VER}/VerifyBamID2.Linux.x86_64.tar.gz"
    if fetch "$VBID_URL" "$TOOL_DIR/src/verifybamid2/vbid.tar.gz" 2>/dev/null; then
        tar -xzf "$TOOL_DIR/src/verifybamid2/vbid.tar.gz" -C "$TOOL_DIR/verifybamid2/" --strip-components=1 || \
        tar -xzf "$TOOL_DIR/src/verifybamid2/vbid.tar.gz" -C "$TOOL_DIR/verifybamid2/"
        # Locate the binary and resource directory
        VBID_BIN=$(find "$TOOL_DIR/verifybamid2" -name "VerifyBamID" -type f | head -1)
        if [[ -n "$VBID_BIN" ]]; then
            chmod +x "$VBID_BIN"
            ln -sf "$VBID_BIN" "$TOOL_DIR/bin/VerifyBamID"
            ok "VerifyBamID2: $TOOL_DIR/bin/VerifyBamID"
            VBID_RESOURCE=$(find "$TOOL_DIR/verifybamid2" -name "*.dat" | head -1 | xargs dirname 2>/dev/null || true)
            if [[ -n "$VBID_RESOURCE" ]]; then
                ok "VerifyBamID2 resource files: $VBID_RESOURCE"
            else
                warn "VerifyBamID2 SVD resource files not found in archive. Download them manually."
                warn "See: https://github.com/Griffan/VerifyBamID/tree/master/resource"
            fi
        fi
    else
        warn "VerifyBamID2 download failed — pipeline will use GATK CalculateContamination fallback."
    fi
fi

# ── 10. Python packages ────────────────────────────────────────────────────────
# The pipeline uses python3 inline scripts. Core packages (sys, math, gzip, os)
# are stdlib. pandas/numpy are only used in final_report.nf.
# If python3 is available, install missing packages to user space.
if command -v python3 &>/dev/null; then
    ok "Python3: $(python3 --version)"
    python3 -c "import pandas" 2>/dev/null || {
        echo "Installing pandas/numpy to user space..."
        python3 -m pip install --user pandas numpy scipy 2>&1 | tail -3
    }
else
    warn "python3 not found. On HPC try:  module load python"
fi

# ── 11. R (optional — used for QC plots) ──────────────────────────────────────
# R and ggplot2 are optional; plots are skipped if Rscript is absent.
if command -v Rscript &>/dev/null; then
    ok "Rscript: $(Rscript --version 2>&1)"
    Rscript -e "library(ggplot2)" 2>/dev/null || {
        echo "Installing ggplot2..."
        Rscript -e "install.packages('ggplot2', repos='https://cloud.r-project.org', lib=Sys.getenv('R_LIBS_USER'))"
    }
else
    warn "Rscript not found — QC plots will be skipped (pipeline still works)."
fi

# ── Summary ────────────────────────────────────────────────────────────────────
echo ""
echo "================================================================"
echo " Setup complete. To use these tools with the pipeline:"
echo "================================================================"
echo ""
echo " Option A — set an environment variable (recommended):"
echo "   export GENETIC_QC_TOOL_DIR=$TOOL_DIR"
echo "   # Then run the pipeline:"
echo "   nextflow run snp_array_qc/main.nf -profile manual_paths ..."
echo "   nextflow run wgs_wes_qc/main.nf   -profile slurm,manual_paths ..."
echo ""
echo " Option B — pass --tool_dir at runtime:"
echo "   nextflow run snp_array_qc/main.nf \\"
echo "     --tool_dir $TOOL_DIR \\"
echo "     -profile manual_paths ..."
echo ""
echo " Verify the tools:"
echo "   export PATH=$TOOL_DIR/bin:\$PATH"
echo "   plink --version && plink2 --version && samtools --version | head -1"
echo "   bcftools --version | head -1 && gatk --version && fastqc --version"
echo "   mosdepth --version"
echo "================================================================"
echo ""
