#!/usr/bin/env bash
set -euo pipefail

IMAGE="${GENETIC_QC_DOCKER_IMAGE:-ghcr.io/kaiyao28/genetic-qc:1.0}"
FORCE_PULL="${GENETIC_QC_FORCE_PULL:-false}"

echo "Genetic QC smoke tests"
echo "Docker image: ${IMAGE}"
echo

if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: docker is not available on PATH."
    echo "Install/start Docker first, then re-run this script."
    exit 1
fi

if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: nextflow is not available on PATH."
    echo "Install Nextflow on the host machine, then re-run this script."
    exit 1
fi

if [ ! -f "nextflow.config" ]; then
    echo "ERROR: run this script from the repository root:"
    echo "  bash test_data/run_smoke_tests.sh"
    exit 1
fi

echo "Docker storage summary:"
docker system df || true
echo

print_docker_recovery_help() {
    cat << EOF

Docker could not pull/extract the image.

This is usually a Docker Desktop / WSL storage problem, not a pipeline problem.
Common causes are:
  - Docker Desktop virtual disk is full
  - a previous failed pull left a partial/broken image layer
  - Docker Desktop needs a restart

Try these steps, then re-run:

  docker image rm ${IMAGE}
  docker builder prune
  docker system prune

If Docker still reports input/output error, restart Docker Desktop.
If it still fails after restart, open Docker Desktop:

  Settings -> Resources -> Advanced -> increase disk image size

or:

  Troubleshoot -> Clean / Purge data

Warning: Docker Desktop purge removes local Docker images and containers, but
not your Git repository files.
EOF
}

echo "Checking Docker image..."
if docker image inspect "${IMAGE}" >/dev/null 2>&1 && [ "${FORCE_PULL}" != "true" ]; then
    echo "Image already exists locally; skipping docker pull."
    echo "Set GENETIC_QC_FORCE_PULL=true to force a fresh pull."
else
    if ! docker pull "${IMAGE}"; then
        print_docker_recovery_help
        exit 1
    fi
fi
echo

echo "Preparing SNP-array PLINK binary test data..."
if [ ! -f "test_data/snp_array/toy.bed" ] || \
   [ ! -f "test_data/snp_array/toy.bim" ] || \
   [ ! -f "test_data/snp_array/toy.fam" ]; then
    docker run --rm \
        -v "${PWD}/test_data/snp_array:/data" \
        "${IMAGE}" \
        bash -lc "cd /data && plink --file toy --make-bed --out toy --allow-no-sex"
else
    echo "SNP-array binary files already exist."
fi
echo

echo "Running WGS/WES VCF variant-only smoke test..."
nextflow run wgs_wes_qc/main.nf \
  --input_type vcf \
  --samplesheet test_data/wgs_wes/samplesheet_vcf.csv \
  --reference_fasta test_data/reference/mini.fa \
  --mode wgs \
  --chroms 22 \
  --run_variant_qc true \
  --run_variant_filtering false \
  --run_sample_qc false \
  --run_final_report true \
  --outdir results/test_vcf_variant_only \
  --docker_image "${IMAGE}" \
  -profile docker \
  -ansi-log false \
  -resume
echo

echo "Running SNP-array variant-only smoke test..."
nextflow run snp_array_qc/main.nf \
  --bfile test_data/snp_array/toy \
  --run_variant_qc true \
  --run_sample_qc false \
  --run_final_report true \
  --outdir results/test_snp_variant_only \
  --docker_image "${IMAGE}" \
  -profile docker \
  -ansi-log false \
  -resume
echo

echo "Smoke tests finished."
echo "Check:"
echo "  results/test_vcf_variant_only"
echo "  results/test_snp_variant_only"
