/*
================================================================================
  MODULE: PCA_VARIANT_SELECTION (WGS/WES only)
================================================================================
  Purpose:
    Select a high-quality, common-variant subset from a sequencing VCF that is
    suitable for PCA. This step is REQUIRED for WGS/WES data and does NOT
    exist in the SNP array pipeline.

  Why WGS/WES needs this step (and SNP arrays do not):
    SNP arrays already genotype common variants by design — after standard QC
    (missingness, HWE) the remaining variants are already suitable for PCA.

    WGS/WES VCFs contain:
      - Rare variants (MAF < 5%) that distort PCA toward rare-variant structure
      - Multi-allelic sites that PLINK cannot use for PCA
      - Indels with unreliable depth estimates
      - Sites with variable coverage (especially WES capture edge effects)
      - Potentially millions of variants — using all is computationally wasteful
        and statistically unsound

    This module produces a clean, common-SNP VCF that feeds directly into
    LD pruning and PCA.

  Additional WES caveat:
    WES PCA may be modestly biased because the capture kit enriches certain
    coding regions. Users should interpret WES PCA components with caution
    and prefer WGS or array data for ancestry assignment when available.
    The pipeline logs a warning when mode = "wes".

  Params controlled here:
    params.pca_maf                 = 0.05   min MAF for PCA variants
    params.pca_variant_missingness = 0.01   max site-level missingness
    params.pca_use_snps_only       = true   exclude indels
    params.pca_autosomes_only      = true   restrict to autosomes (recommended)

  How to change:
    nextflow run main.nf --pca_maf 0.01    # include less-common variants
    nextflow run main.nf --pca_autosomes_only false  # include chrX

  Output:
    - pca_variants.vcf.gz       : filtered VCF for PCA (common biallelic SNPs)
    - pca_variants.vcf.gz.tbi   : tabix index
    - pca_variant_selection_summary.txt : counts before/after filtering
================================================================================
*/

process PCA_VARIANT_SELECTION {
    label 'process_medium'
    publishDir "${params.outdir}/ancestry_pca", mode: params.publish_dir_mode, pattern: "pca_variant_selection_summary.txt"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("pca_variants.vcf.gz"),
                     path("pca_variants.vcf.gz.tbi"), emit: vcf
    path "pca_variant_selection_summary.txt",          emit: summary

    script:
    def snps_only_flag  = params.pca_use_snps_only   ? "--type snps"           : ""
    def auto_only_flag  = params.pca_autosomes_only   ? "--regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22" : ""
    """
    # ── Warn if WES mode (PCA may be biased) ─────────────────────────────────
    if [ "${meta.mode}" = "wes" ]; then
        echo "WARNING: WES PCA variant selection is restricted to capture regions."
        echo "         Coverage biases from probe efficiency may affect PC interpretation."
        echo "         Prefer WGS or SNP array data for ancestry assignment when available."
    fi

    n_before=\$(bcftools view --no-header ${vcf} | wc -l)

    # ── Step 1: biallelic SNPs only ───────────────────────────────────────────
    # Multi-allelic sites and indels are unreliable for PCA
    bcftools view \\
        --threads ${task.cpus} \\
        ${snps_only_flag} \\
        --min-alleles 2 --max-alleles 2 \\
        ${auto_only_flag} \\
        ${vcf} \\
        -O z -o step1_biallelic_snps.vcf.gz

    bcftools index --tbi step1_biallelic_snps.vcf.gz

    # ── Step 2: PASS filter ───────────────────────────────────────────────────
    # Only use variants that passed site-level quality filters
    bcftools view \\
        --threads ${task.cpus} \\
        --apply-filters PASS \\
        step1_biallelic_snps.vcf.gz \\
        -O z -o step2_pass.vcf.gz

    bcftools index --tbi step2_pass.vcf.gz

    # ── Step 3: MAF ≥ pca_maf ────────────────────────────────────────────────
    # Rare variants distort PCA toward rare-variant structure and are poorly
    # estimated — exclude everything below the PCA MAF threshold.
    bcftools view \\
        --threads ${task.cpus} \\
        --min-af ${params.pca_maf}:minor \\
        step2_pass.vcf.gz \\
        -O z -o step3_maf.vcf.gz

    bcftools index --tbi step3_maf.vcf.gz

    # ── Step 4: low missingness ───────────────────────────────────────────────
    # Sites with missing data in many samples produce unreliable PC scores
    bcftools filter \\
        --threads ${task.cpus} \\
        --include "F_MISSING < ${params.pca_variant_missingness}" \\
        step3_maf.vcf.gz \\
        -O z -o pca_variants.vcf.gz

    bcftools index --tbi pca_variants.vcf.gz

    n_after=\$(bcftools view --no-header pca_variants.vcf.gz | wc -l)
    n_removed=\$(( n_before - n_after ))

    # ── Summary ───────────────────────────────────────────────────────────────
    cat > pca_variant_selection_summary.txt << EOF
step=pca_variant_selection
dataset=${meta.id}
mode=${meta.mode}
pca_maf=${params.pca_maf}
pca_variant_missingness=${params.pca_variant_missingness}
pca_use_snps_only=${params.pca_use_snps_only}
pca_autosomes_only=${params.pca_autosomes_only}
n_variants_before=\${n_before}
n_variants_removed=\${n_removed}
n_variants_selected=\${n_after}
EOF

    echo "PCA variant selection: \${n_after} variants selected from \${n_before} (removed \${n_removed})"
    echo "  Filters: biallelic SNPs, PASS, MAF >= ${params.pca_maf}, F_MISSING < ${params.pca_variant_missingness}"
    """
}
