# Test Data

Small toy files for checking that the pipeline starts, parses inputs, and reaches the expected modules.

These files are intentionally tiny. They are for smoke testing pipeline wiring, not for biological interpretation.

## Run All Smoke Tests

From the repository root, run:

```bash
bash test_data/run_smoke_tests.sh
```

This script checks that Docker and Nextflow are available, pulls the published Docker image, creates the SNP-array PLINK binary test files inside Docker, then runs:

```text
WGS/WES VCF variant-only smoke test
SNP-array variant-only smoke test
```

It is the recommended first check after cloning the repository.

If the Docker image is already present locally, the script skips `docker pull`.
To force a fresh pull:

```bash
GENETIC_QC_FORCE_PULL=true bash test_data/run_smoke_tests.sh
```

If Docker fails with an `input/output error` while extracting an image layer,
restart Docker Desktop and clean old Docker data:

```bash
docker image rm ghcr.io/kaiyao28/genetic-qc:1.0
docker builder prune
docker system prune
```

## WGS/WES VCF Smoke Test

Use the VCF fixture for a quick variant-level test:

```bash
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
  -profile docker
```

This should exercise:

```text
01_input_check
02_variant_level_qc
04_final_report
```

For a sample-level smoke test on this tiny VCF, keep expectations modest and disable modules that need many independent variants:

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type vcf \
  --samplesheet test_data/wgs_wes/samplesheet_vcf.csv \
  --reference_fasta test_data/reference/mini.fa \
  --mode wgs \
  --chroms 22 \
  --sample_qc_scope provisional \
  --run_relatedness_wgs false \
  --run_ancestry_pca_wgs false \
  --run_final_report true \
  --outdir results/test_vcf_sample_smoke \
  -profile docker
```

## SNP Array Smoke Test

The SNP-array fixture is provided as PED/MAP text files so the repository does not need to store binary PLINK files. Convert it before running the SNP-array pipeline:

```bash
docker run --rm \
  -v "$PWD/test_data/snp_array:/data" \
  ghcr.io/kaiyao28/genetic-qc:1.0 \
  bash -lc "cd /data && plink --file toy --make-bed --out toy --allow-no-sex"
```

Then run:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile test_data/snp_array/toy \
  --run_variant_qc true \
  --run_sample_qc false \
  --run_final_report true \
  --outdir results/test_snp_variant_only \
  -profile docker
```

The Docker conversion avoids needing `plink` on the host machine. If `plink` is already installed locally, `bash test_data/snp_array/make_plink_binary.sh` also works.
