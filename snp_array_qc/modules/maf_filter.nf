/*
================================================================================
  MODULE: MAF_FILTER
================================================================================
  Purpose:
    Remove variants with minor allele frequency (MAF) below a threshold.
    Rare variants have insufficient power for standard GWAS-style analyses and
    their call rates and HWE tests are less reliable.

  Why this step:
    For common-variant GWAS and polygenic score applications, MAF ≥ 0.01 (1%)
    is standard. Very rare variants (MAF < 0.001) are poorly estimated on SNP
    arrays and should almost always be removed in array-based studies.

    For rare-variant association tests (burden, SKAT), set:
      params.run_maf_filter = false
    or use a much lower threshold (e.g. params.maf = 0.001).

  Default threshold: params.maf = 0.01

  How to change:
    nextflow run main.nf --maf 0.05    # more stringent (common variants only)
    nextflow run main.nf --maf 0.001   # keep rare variants

  How to disable:
    nextflow run main.nf --run_maf_filter false

  Output:
    - Filtered PLINK files
    - maf_removed_variants.txt   : IDs of removed variants
    - maf_summary.txt            : counts before/after
================================================================================
*/

process MAF_FILTER {
    label 'process_medium'
    publishDir "${params.outdir}/cleaned_data", mode: params.publish_dir_mode, pattern: "*.{bed,bim,fam}"
    publishDir "${params.outdir}/qc_tables",    mode: params.publish_dir_mode, pattern: "*.txt"
    publishDir "${params.outdir}/qc_plots",    mode: params.publish_dir_mode, pattern: "*.png"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("${meta.id}_maf.bed"),
                     path("${meta.id}_maf.bim"),
                     path("${meta.id}_maf.fam"), emit: plink
    path "maf_removed_variants.txt",              emit: removed_variants
    path "maf_summary.txt",                       emit: summary
    path "*.png",                                 optional: true, emit: plots

    script:
    def prefix = "${meta.id}"
    """
    n_var_before=\$(wc -l < ${bim})

    plink \\
        --bfile ${bed.baseName} \\
        --maf ${params.maf} \\
        --make-bed \\
        --out ${prefix}_maf \\
        --allow-no-sex

    n_var_after=\$(wc -l < ${prefix}_maf.bim)
    n_removed=\$(( n_var_before - n_var_after ))

    awk 'NR==FNR{seen[\$2]=1; next} !seen[\$2]{print \$2}' \\
        ${prefix}_maf.bim ${bim} > maf_removed_variants.txt

    cat > maf_summary.txt << EOF
step=maf_filter
dataset=${meta.id}
maf_threshold=${params.maf}
n_variants_before=\${n_var_before}
n_variants_removed=\${n_removed}
n_variants_after=\${n_var_after}
EOF

    echo "MAF filter (MAF >= ${params.maf}): removed \${n_removed} variants"

    # ── MAF distribution plot (computed on pre-filter data) ───────────────────
    plink \\
        --bfile ${bed.baseName} \\
        --freq \\
        --out maf_dist \\
        --allow-no-sex

    if command -v Rscript &>/dev/null; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
df <- read.table("maf_dist.frq", header=TRUE)
p <- ggplot(df, aes(x=MAF)) +
    geom_histogram(bins=100, fill="steelblue", alpha=0.8) +
    geom_vline(xintercept=${params.maf}, linetype="dashed", colour="red") +
    labs(title="Minor allele frequency (MAF) distribution",
         subtitle="Red line: removal threshold — variants to the left are removed",
         x="Minor allele frequency", y="Number of variants") +
    theme_classic()
ggsave("maf_plot.png", p, width=8, height=5)
RSCRIPT
    fi
    """
}
