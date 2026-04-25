# Test Data

Small toy files for checking that the pipeline starts, parses inputs, and reaches the expected modules.

These files are intentionally tiny. They are for smoke testing pipeline wiring, not for biological interpretation.

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
  --run_sample_qc false \
  --run_final_report true \
  --outdir results/test_vcf_variant_only \
  -profile standard
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
  -profile standard
```

## SNP Array Smoke Test

The SNP-array fixture is provided as PED/MAP text files so the repository does not need to store binary PLINK files. Convert it before running the SNP-array pipeline:

```bash
cd test_data/snp_array
bash make_plink_binary.sh
cd ../..
```

Then run:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile test_data/snp_array/toy \
  --run_variant_qc true \
  --run_sample_qc false \
  --run_final_report true \
  --outdir results/test_snp_variant_only \
  -profile standard
```

The conversion step requires `plink` on `PATH`.
