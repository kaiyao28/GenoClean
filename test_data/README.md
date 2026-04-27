# Test Data

Small synthetic files for smoke-testing pipeline wiring and QC filter logic.
These files are intentionally tiny — biological interpretation of results is meaningless.

## Regenerating Test Data

If the test data files are missing or you want to reset them:

```bash
python3 test_data/generate_test_data.py
rm -f test_data/snp_array/toy.bed test_data/snp_array/toy.bim test_data/snp_array/toy.fam
```

## Running Smoke Tests

Run both tests (default):

```bash
bash test_data/run_smoke_tests.sh --profile manual_paths
bash test_data/run_smoke_tests.sh --profile docker
bash test_data/run_smoke_tests.sh --profile slurm,singularity
```

Run only one test at a time:

```bash
bash test_data/run_smoke_tests.sh --profile manual_paths --test snp_array
bash test_data/run_smoke_tests.sh --profile manual_paths --test wgs_wes
```

## What Each Test Exercises

### SNP Array (`--test snp_array`)

Input: `test_data/snp_array/toy.ped` and `toy.map`

```text
10 samples (5 male, 5 female, 6 cases, 4 controls)
50 variants on chr22

Variants 1–40  : normal genotypes, MAF 0.10–0.40, no missingness
Variants 41–45 : ~35% missing genotypes → fail variant callrate filter (>2%)
Variants 46–48 : monomorphic (MAF=0) → fail MAF filter (<0.01)
Variants 49–50 : all samples heterozygous → exercises HWE step
```

Run command used in the smoke test:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile test_data/snp_array/toy \
  --run_variant_qc true \
  --run_sample_qc false \
  --run_final_report true \
  --outdir results/test_snp_variant_only \
  -profile manual_paths
```

### WGS/WES (`--test wgs_wes`)

Input: `test_data/wgs_wes/toy_chr22.vcf` (VCF mode, variant QC only, no filtering)

```text
5 samples (SAMPLE1–5)
20 variants on chr22:
  18 SNPs: 13 transitions + 5 transversions → Ti/Tv ≈ 2.6
  2 indels (1 deletion, 1 insertion)
  rs013 (pos 325): FAIL — QD < 2.0 (low quality by depth)
  rs015 (pos 375): FAIL — FS > 60.0 (strand bias)
  rs018 (pos 450): PASS site, but 2 samples have GQ < 20 (genotype filter demo)
```

Reference: `test_data/reference/mini.fa` — 500bp deterministic chr22 (ATCGATCG repeating).
REF alleles in the VCF are verified to match this reference at each position.

Run command used in the smoke test:

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
  -profile manual_paths
```
