# Genetic QC Nextflow Pipeline

A reproducible Nextflow DSL2 pipeline for genetic quality control.

The project contains two independent workflows:

- `snp_array_qc/`: QC for PLINK binary files from SNP-array or GWAS-style data.
- `wgs_wes_qc/`: QC for whole-genome or whole-exome sequencing data.

The pipeline separates variant-level QC from sample-level QC. This is important because variant-level QC can often be run per chromosome, while sample-level QC usually requires genome-wide context.

All thresholds are defined in config files and can be overridden at runtime. Individual modules and whole QC phases can be enabled or disabled with parameters.

## Quick Start

Clone the repository and enter the project folder:

```bash
git clone https://github.com/kaiyao28/GeneticQC.git
cd GeneticQC
```

Recommended setup: use Nextflow plus the published Docker image. Docker provides the bioinformatics tools, but Nextflow still needs to be installed on the host machine.

```bash
docker pull ghcr.io/kaiyao28/genetic-qc:1.0
```

Then run the bundled smoke tests with the Docker profile:

```bash
bash test_data/run_smoke_tests.sh
```

For detailed setup instructions on Windows, Linux/macOS, and HPC clusters, see [Setup Guide](docs/setup.md).

## Workflow Design

Both workflows follow the same high-level structure:

```text
01_input_check
02_variant_level_qc
03_sample_level_qc
04_final_report
```

`04_final_report` is an aggregation step, not another QC filter. It can run after a full QC workflow, after variant-level QC only, or after sample-level QC only. The report records which phases ran, which phases were skipped, which chromosomes were analysed, and whether sample-level conclusions are final or provisional.

Recommended phase combinations:

| Use case | Parameters | Expected report status |
|----------|------------|------------------------|
| Full production QC | `--run_variant_qc true --run_sample_qc true --chroms 1-22 --sample_qc_scope auto` | Variant QC complete; sample QC final; suitable for final filtering. |
| Variant-only QC | `--run_variant_qc true --run_sample_qc false --run_final_report true` | Variant QC complete; sample QC skipped; report still generated. |
| Chromosome test run | `--chroms 22 --sample_qc_scope provisional --run_final_report true` | Variant QC for chr22; sample QC marked provisional if run. |
| Skip all sample QC | `--sample_qc_scope skip --run_final_report true` | Sample-level QC skipped by scope; report still generated. |
| Report disabled | `--run_final_report false` | QC runs, but no final HTML report is produced. |

Variant-level QC includes filters and metrics applied to variants, sites, or genotypes. These steps can often run per chromosome:

- duplicate variant checks
- variant missingness
- HWE filtering
- MAF filtering
- site-level filtering
- genotype-level filtering
- Ti/Tv and variant count summaries
- chromosome merge or aggregation

Sample-level QC includes filters and metrics applied to individuals. These steps should be treated as final only when all autosomes are present:

- sample missingness
- heterozygosity
- sex check
- relatedness
- ancestry PCA
- contamination
- sample-level variant counts
- coverage and alignment metrics

The `sample_qc_scope` parameter controls this behavior:

| Value | Behavior |
|-------|----------|
| `auto` | Use `genome_wide` when all autosomes are present; otherwise mark sample QC as provisional. |
| `genome_wide` | Treat sample-level QC as final. |
| `provisional` | Run sample-level QC but mark it as not suitable for final filtering. |
| `skip` | Skip sample-level QC. |

## SNP Array QC

This example runs the full SNP-array workflow:

```text
01_input_check
02_variant_level_qc
03_sample_level_qc
04_final_report
```

The reference panel is optional. If supplied, it should be a PLINK binary prefix with matching `.bed`, `.bim`, and `.fam` files, for example `data/1000G/1000G_hg38`.

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --reference_panel data/1000G/1000G_hg38 \
  --run_variant_qc true \
  --run_sample_qc true \
  --chroms 1-22 \
  --sample_qc_scope auto \
  --run_final_report true \
  --outdir results/snp_array_qc \
  -profile docker
```

For chromosome-only tests, variant-only runs, threshold overrides, skipped modules, and alternative reference-panel setups, see the [SNP Array QC Manual](docs/snp_array_qc_manual.md).

## WGS / WES QC

This example runs the full WES workflow from BAM/CRAM input:

```text
01_input_check
02_variant_level_qc
03_sample_level_qc
04_final_report
```

BAM/CRAM input supports sample-level sequencing QC such as alignment metrics, duplicate rate, coverage, contamination, and sex check. VCF input supports variant-level QC and VCF-based sample QC such as sample call rate, relatedness, and ancestry PCA.

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type bam \
  --samplesheet assets/example_samplesheet.csv \
  --reference_fasta /data/reference/GRCh38.fa \
  --target_intervals /data/reference/exome_targets.bed \
  --mode wes \
  --run_variant_qc true \
  --run_sample_qc true \
  --chroms 1-22 \
  --sample_qc_scope auto \
  --run_final_report true \
  --outdir results/wgs_wes_qc \
  -profile singularity
```

This example runs the full VCF-based WGS/WES QC path: chromosome-level variant QC, merge, genome-wide sample QC, then final report.

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type vcf \
  --samplesheet samplesheet.csv \
  --reference_fasta /data/reference/GRCh38.fa \
  --mode wgs \
  --chroms 1-22 \
  --run_variant_qc true \
  --run_sample_qc true \
  --sample_qc_scope auto \
  --run_final_report true \
  --outdir results/vcf_qc \
  -profile docker
```

For FASTQ input, VCF-only runs, chromosome-only tests, skipped modules, provisional sample QC, and threshold overrides, see the [WGS/WES QC Manual](docs/wgs_wes_qc_manual.md).

## Params Files

Run with a params file:

```bash
nextflow run snp_array_qc/main.nf -params-file assets/example_params.yml
```

See [assets/example_params.yml](assets/example_params.yml) for an annotated example.

## Changing Defaults

Use command-line flags for one-off changes:

```bash
nextflow run snp_array_qc/main.nf --bfile data/raw/genotypes --maf 0.05
```

Use workflow-specific config files for project defaults:

```text
conf/snp_array_qc.config
conf/wgs_wes_qc.config
```

Use the root `nextflow.config` for shared defaults and execution profiles such as Docker, Singularity, Conda, SLURM, LSF, and AWS Batch.

In general:

```text
temporary run change     -> command-line flags
SNP-array defaults       -> conf/snp_array_qc.config
WGS/WES defaults         -> conf/wgs_wes_qc.config
shared profiles/defaults -> nextflow.config
```

## Test Data

Small toy files for smoke testing are available in [test_data/](test_data/). These examples are intentionally tiny and are meant to check pipeline wiring, not biological validity.

Recommended all-in-one check after cloning:

```bash
bash test_data/run_smoke_tests.sh
```

This checks Docker and Nextflow, pulls `ghcr.io/kaiyao28/genetic-qc:1.0`, prepares the SNP-array example data inside Docker, then runs both example workflows in variant-only report mode.

The script shows Docker disk usage first and skips `docker pull` if the image is already present. To force a fresh pull:

```bash
GENETIC_QC_FORCE_PULL=true bash test_data/run_smoke_tests.sh
```

Quick VCF smoke test:

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

Quick SNP-array smoke test:

```bash
docker run --rm \
  -v "$PWD/test_data/snp_array:/data" \
  ghcr.io/kaiyao28/genetic-qc:1.0 \
  bash -lc "cd /data && plink --file toy --make-bed --out toy --allow-no-sex"

nextflow run snp_array_qc/main.nf \
  --bfile test_data/snp_array/toy \
  --run_variant_qc true \
  --run_sample_qc false \
  --run_final_report true \
  --outdir results/test_snp_variant_only \
  -profile docker
```

See [test_data/README.md](test_data/README.md) for notes and additional smoke-test commands.

## Execution Profiles

Profiles can be combined, for example `-profile slurm,singularity`.

| Profile | When to use |
|---------|-------------|
| `standard` | Tools are already available on `PATH`. |
| `docker` | Local machine or cloud execution with Docker. |
| `singularity` | HPC clusters using Singularity or Apptainer. |
| `conda` | Systems without a container daemon. |
| `slurm` | SLURM scheduler; usually combine with a container profile. |
| `lsf` | IBM LSF scheduler. |
| `awsbatch` | AWS Batch. |

Edit queue and account settings in `nextflow.config`:

```groovy
slurm {
    process.queue          = 'your_queue'
    process.clusterOptions = '--account=your_account'
}
```

## Key Parameters

### Shared

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_variant_qc` | true | Enable variant-level QC phase. |
| `--run_sample_qc` | true | Enable sample-level QC phase. |
| `--chroms` | 1-22 | Chromosomes to process. Accepts `1-22`, `22`, `1,2,22`, or `all`. |
| `--sample_qc_scope` | auto | One of `auto`, `genome_wide`, `provisional`, or `skip`. |
| `--keep_intermediate` | false | Retain intermediate files. |
| `--outdir` | results | Output directory. |

### SNP Array

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bfile` | required | PLINK binary prefix without `.bed`, `.bim`, `.fam`. |
| `--pheno` | null | Optional phenotype or covariate file. |
| `--reference_panel` | null | Optional reference panel prefix for ancestry PCA. |
| `--sample_missingness` | 0.02 | Maximum missing genotype fraction per sample. |
| `--variant_missingness` | 0.02 | Maximum missing genotype fraction per variant. |
| `--hwe_p` | 1e-6 | HWE p-value cutoff. |
| `--maf` | 0.01 | Minimum minor allele frequency. |
| `--heterozygosity_sd` | 3 | Heterozygosity outlier SD threshold. |
| `--relatedness_pi_hat` | 0.1875 | Relatedness cutoff. |
| `--pca_outlier_sd` | 6 | PCA outlier SD threshold. |

### WGS / WES

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input_type` | bam | One of `fastq`, `bam`, `cram`, or `vcf`. |
| `--mode` | wes | One of `wes` or `wgs`. |
| `--samplesheet` | required | CSV with `sample,file1[,file2]`. |
| `--reference_fasta` | required | Reference FASTA. A `.fai` index is expected. |
| `--target_intervals` | required for WES | Capture target BED file. |
| `--min_mean_depth_wes` | 30 | Minimum mean on-target depth. |
| `--min_mean_depth_wgs` | 20 | Minimum mean genome-wide depth. |
| `--min_target_20x_fraction` | 0.80 | Minimum fraction of targets covered at least 20x. |
| `--max_contamination` | 0.03 | Maximum allowed contamination estimate. |
| `--max_duplication_rate` | 0.20 | Maximum allowed duplicate rate. |
| `--min_gq` | 20 | Minimum genotype quality. |
| `--min_dp` | 10 | Minimum genotype depth. |
| `--variant_filter_method` | hard_filter | One of `hard_filter` or `vqsr`. |
| `--pca_maf` | 0.05 | Minimum MAF for PCA variant selection. |
| `--pca_autosomes_only` | true | Restrict PCA variant selection to autosomes. |

## Output Structure

SNP-array output:

```text
results/snp_array_qc/
|-- logs/
|-- cleaned_data/
|-- exclusion_lists/
|-- qc_tables/
|-- qc_plots/
`-- final_report.html
```

WGS/WES output:

```text
results/wgs_wes_qc/
|-- logs/
|-- fastqc/
|-- alignment_metrics/
|-- duplicate_metrics/
|-- coverage_qc/
|-- contamination/
|-- sex_check/
|-- variant_calling_qc/
|-- ancestry_pca/
|-- cleaned_data/
|-- qc_tables/
|-- qc_plots/
`-- wgs_wes_final_report.html
```

## Documentation

- [SNP Array QC Manual](docs/snp_array_qc_manual.md)
- [WGS/WES QC Manual](docs/wgs_wes_qc_manual.md)
- [References](docs/references.md)

## References

- Anderson et al. 2010, *Nature Protocols*: GWAS QC thresholds.
- GATK Best Practices: WGS/WES variant filtering and VQSR.
- gnomAD v4.0: contamination cutoffs and sample-level QC framework.
- PLINK 1.9 and PLINK 2 documentation.
