#!/usr/bin/env bash
# =============================================================================
#  test_env.sh — verify all pipeline tools are installed and callable
# =============================================================================
#  Usage:
#    bash test_env.sh            # test tools on PATH (conda / standard mode)
#    bash test_env.sh docker     # test tools inside the Docker image
#    bash test_env.sh singularity # test tools inside the SIF
#
#  Output:
#    Prints PASS / WARN / FAIL per tool to the terminal.
#    Writes a full log to:  test_results.log
#    Exit code 0 = all required tools pass.
#    Exit code 1 = at least one required tool is missing.
# =============================================================================
set -uo pipefail

MODE="${1:-path}"
LOG="test_results.log"
IMAGE="${GENETIC_QC_DOCKER_IMAGE:-ghcr.io/kaiyao28/genetic-qc:1.0}"
SIF="containers/genetic-qc.sif"

# ── colours ───────────────────────────────────────────────────────────────────
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'

# ── run a single command, optionally inside a container ───────────────────────
run_cmd() {
    if [[ "$MODE" == "docker" ]]; then
        docker run --rm "$IMAGE" bash -c "$*" 2>&1
    elif [[ "$MODE" == "singularity" ]]; then
        apptainer exec "$SIF" bash -c "$*" 2>&1
    else
        bash -c "$*" 2>&1
    fi
}

# ── check one tool ────────────────────────────────────────────────────────────
#  check_tool LABEL REQUIRED VERSION_CMD
PASS_COUNT=0
WARN_COUNT=0
FAIL_COUNT=0

check_tool() {
    local label="$1"
    local required="$2"   # "required" or "optional"
    local cmd="$3"

    local out
    out=$(run_cmd "$cmd" 2>&1 | head -1) && local rc=0 || local rc=$?

    if [[ -n "$out" && "$rc" -eq 0 ]] || [[ -n "$out" && "$required" == "optional" ]]; then
        printf "${GREEN}[PASS]${NC}  %-28s %s\n" "$label" "$out"
        echo "[PASS]  ${label}: ${out}" >> "$LOG"
        (( PASS_COUNT++ )) || true
    elif [[ "$required" == "optional" ]]; then
        printf "${YELLOW}[WARN]${NC}  %-28s not found (optional — some features disabled)\n" "$label"
        echo "[WARN]  ${label}: not found (optional)" >> "$LOG"
        (( WARN_COUNT++ )) || true
    else
        printf "${RED}[FAIL]${NC}  %-28s NOT FOUND — required\n" "$label"
        echo "[FAIL]  ${label}: NOT FOUND — required" >> "$LOG"
        (( FAIL_COUNT++ )) || true
    fi
}

# ── header ────────────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}=== Genetic QC Pipeline — Environment Test (mode: ${MODE}) ===${NC}"
echo "$(date)" | tee "$LOG"
echo "mode: ${MODE}" >> "$LOG"
echo "" | tee -a "$LOG"

# ── Nextflow ──────────────────────────────────────────────────────────────────
NF_VER=$(nextflow -version 2>&1 | grep -oP 'version \K[\d.]+' | head -1 || echo "")
if [[ -n "$NF_VER" ]]; then
    printf "${GREEN}[PASS]${NC}  %-28s nextflow %s\n" "Nextflow" "$NF_VER"
    echo "[PASS]  Nextflow: ${NF_VER}" >> "$LOG"
    (( PASS_COUNT++ )) || true
else
    printf "${RED}[FAIL]${NC}  %-28s NOT FOUND — required\n" "Nextflow"
    echo "[FAIL]  Nextflow: NOT FOUND" >> "$LOG"
    (( FAIL_COUNT++ )) || true
fi

# ── Java (required by Nextflow, GATK, Picard) ─────────────────────────────────
check_tool "Java"         required  "java -version 2>&1 | head -1"

# ── SNP array tools ───────────────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo -e "${BOLD}--- SNP array tools ---${NC}"
echo "--- SNP array tools ---" >> "$LOG"

check_tool "PLINK 1.9"    required  "plink --version 2>&1 | head -1"
check_tool "PLINK2"       required  "plink2 --version 2>&1 | head -1"

# ── WGS/WES tools ─────────────────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo -e "${BOLD}--- WGS/WES tools ---${NC}"
echo "--- WGS/WES tools ---" >> "$LOG"

check_tool "samtools"     required  "samtools --version 2>&1 | head -1"
check_tool "bcftools"     required  "bcftools --version 2>&1 | head -1"
check_tool "picard"       required  "picard --version 2>&1 | head -1"
check_tool "FastQC"       required  "fastqc --version 2>&1 | head -1"
check_tool "mosdepth"     required  "mosdepth --version 2>&1 | head -1"
check_tool "GATK"         required  "gatk --version 2>&1 | head -1"
check_tool "tabix"        required  "tabix --version 2>&1 | head -1"
check_tool "VerifyBamID2" optional  "VerifyBamID --version 2>&1 | head -1"

# ── Python + key packages ─────────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo -e "${BOLD}--- Python ---${NC}"
echo "--- Python ---" >> "$LOG"

check_tool "Python 3"     required  "python3 --version 2>&1"
check_tool "  pandas"     required  "python3 -c 'import pandas; print(pandas.__version__)'"
check_tool "  numpy"      required  "python3 -c 'import numpy;  print(numpy.__version__)'"
check_tool "  scipy"      required  "python3 -c 'import scipy;  print(scipy.__version__)'"

# ── R + key packages (optional) ───────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo -e "${BOLD}--- R (optional — used for plots) ---${NC}"
echo "--- R ---" >> "$LOG"

check_tool "R"            optional  "Rscript --version 2>&1 | head -1"
check_tool "  ggplot2"    optional  "Rscript -e 'cat(as.character(packageVersion(\"ggplot2\")))'"

# ── Summary ───────────────────────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo "================================================" | tee -a "$LOG"
printf "${BOLD}Results:${NC}  "
printf "${GREEN}%d PASS${NC}  " "$PASS_COUNT"
[[ $WARN_COUNT -gt 0 ]] && printf "${YELLOW}%d WARN${NC}  " "$WARN_COUNT"
[[ $FAIL_COUNT -gt 0 ]] && printf "${RED}%d FAIL${NC}"    "$FAIL_COUNT"
echo ""
echo "Results: ${PASS_COUNT} PASS  ${WARN_COUNT} WARN  ${FAIL_COUNT} FAIL" >> "$LOG"
echo "Log written to: ${LOG}"

if [[ $FAIL_COUNT -gt 0 ]]; then
    echo ""
    echo -e "${RED}FAIL: ${FAIL_COUNT} required tool(s) missing.${NC}"
    echo "Check ${LOG} for details, then re-run setup.sh."
    echo ""
    exit 1
elif [[ $WARN_COUNT -gt 0 ]]; then
    echo ""
    echo -e "${YELLOW}WARN: optional tools missing — pipeline will run but plots/contamination check may be skipped.${NC}"
    echo ""
    exit 0
else
    echo ""
    echo -e "${GREEN}All checks passed. Ready to run the pipeline.${NC}"
    echo ""
    exit 0
fi
