/*
================================================================================
  MODULE: HWE_FILTER
================================================================================
  Purpose:
    Remove SNPs that significantly deviate from Hardy-Weinberg equilibrium (HWE).
    Extreme HWE deviation in a general population usually indicates genotyping
    error rather than true biology.

  Why this step:
    HWE deviation can arise from:
      - Genotyping errors (most common cause)
      - True selection at the locus (rare for common variants)
      - Population stratification (use controls or account for ancestry first)
      - Null alleles or probe-binding failures on arrays

    IMPORTANT: In case-control studies, HWE should be evaluated in CONTROLS ONLY.
    HWE deviation in cases may reflect genuine disease association, not error.
    This module applies HWE filtering in controls if phenotype codes are present;
    otherwise it applies the filter across all samples.

  Default threshold: params.hwe_p = 1e-6

  How to change:
    nextflow run main.nf --hwe_p 1e-4   # more permissive
    nextflow run main.nf --hwe_p 1e-10  # more stringent

  How to disable:
    nextflow run main.nf --run_hwe false

  Output:
    - Filtered PLINK files
    - hwe_removed_variants.txt   : IDs of removed variants
    - hwe_summary.txt            : counts before/after
================================================================================
*/

process HWE_FILTER {
    label 'process_medium'
    publishDir "${params.outdir}/qc_tables",   mode: params.publish_dir_mode, pattern: "*.txt"
    publishDir "${params.outdir}/qc_plots",   mode: params.publish_dir_mode, pattern: "*.png"
    publishDir "${params.outdir}/cleaned_data", mode: params.publish_dir_mode, pattern: "*.{bed,bim,fam}", enabled: params.keep_intermediate

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path pheno   // may be empty; used to detect case-control labels

    output:
    tuple val(meta), path("${meta.id}_hwe.bed"),
                     path("${meta.id}_hwe.bim"),
                     path("${meta.id}_hwe.fam"), emit: plink
    path "hwe_removed_variants.txt",              emit: removed_variants
    path "hwe_summary.txt",                       emit: summary
    path "*.png",                                 optional: true, emit: plots

    script:
    def prefix = "${meta.id}"
    """
    n_var_before=\$(wc -l < ${bim})

    # ── Detect case-control phenotype coding ──────────────────────────────────
    # FAM column 6: 1=control, 2=case; 0/-9=missing
    n_cases=\$(awk '\$6==2' ${fam} | wc -l)
    n_controls=\$(awk '\$6==1' ${fam} | wc -l)

    if [ "\${n_controls}" -gt 0 ]; then
        echo "Case-control data: \${n_cases} cases, \${n_controls} controls — HWE in CONTROLS ONLY"
        hwe_test="UNAFF"
        filter_mode="controls_only"
    else
        echo "No case-control coding — HWE across all samples"
        hwe_test="ALL"
        filter_mode="all_samples"
    fi

    # ── Compute HWE statistics on all variants ────────────────────────────────
    # --hardy produces a .hwe file with TEST column: ALL / AFF / UNAFF
    # Using this gives us per-chromosome control over thresholds, which
    # plink --hwe does not support (it applies one threshold to all chromosomes)
    plink \\
        --bfile ${bed.baseName} \\
        --hardy \\
        --out hwe_stats \\
        --allow-no-sex

    # ── Build variant exclusion list with chromosome-aware thresholds ─────────
    # Detect whether chrX (CHR 23) SNPs are present in this dataset.
    # Arrays restricted to autosomes (1-22) are common; in that case the
    # separate chrX threshold is not applicable and the autosome threshold
    # is used for all variants.
    n_chrx_snps=\$(awk 'NR>1 && \$1+0==23' hwe_stats.hwe | wc -l)

    if [ "\${n_chrx_snps}" -gt 0 ]; then
        # chrX present — apply separate, more permissive threshold.
        # Males are hemizygous on X; PLINK computes X HWE in females only,
        # making the test noisier. Using a relaxed threshold avoids over-filtering.
        awk -v a="${params.hwe_p}" -v x="${params.hwe_p_chrx}" -v t="\${hwe_test}" \\
            'NR>1 && \$3==t && ((\$1+0<23 && \$9+0<a) || (\$1+0==23 && \$9+0<x)) {print \$2}' \\
            hwe_stats.hwe > hwe_exclude.txt
        chrx_mode="separate_threshold_p${params.hwe_p_chrx}"
        echo "chrX SNPs detected (\${n_chrx_snps}) — using chrX threshold p<${params.hwe_p_chrx}"
    else
        # No chrX SNPs — autosome threshold applied to all variants
        awk -v a="${params.hwe_p}" -v t="\${hwe_test}" \\
            'NR>1 && \$3==t && \$9+0<a {print \$2}' \\
            hwe_stats.hwe > hwe_exclude.txt
        chrx_mode="no_chrx_snps"
        echo "No chrX SNPs in dataset — autosome threshold applied to all variants"
    fi

    n_hwe_fail=\$(wc -l < hwe_exclude.txt)

    # ── Apply filter ──────────────────────────────────────────────────────────
    if [ "\${n_hwe_fail}" -gt 0 ]; then
        plink \\
            --bfile ${bed.baseName} \\
            --exclude hwe_exclude.txt \\
            --make-bed \\
            --out ${prefix}_hwe \\
            --allow-no-sex
    else
        plink \\
            --bfile ${bed.baseName} \\
            --make-bed \\
            --out ${prefix}_hwe \\
            --allow-no-sex
    fi

    n_var_after=\$(wc -l < ${prefix}_hwe.bim)
    n_removed=\$(( n_var_before - n_var_after ))

    # ── Record removed variant IDs ────────────────────────────────────────────
    awk 'NR==FNR{seen[\$2]=1; next} !seen[\$2]{print \$2}' \\
        ${prefix}_hwe.bim ${bim} > hwe_removed_variants.txt

    cat > hwe_summary.txt << EOF
step=hwe_filter
dataset=${meta.id}
hwe_p_autosomes=${params.hwe_p}
hwe_p_chrx=${params.hwe_p_chrx}
chrx_mode=\${chrx_mode}
filter_mode=\${filter_mode}
n_variants_before=\${n_var_before}
n_variants_removed=\${n_removed}
n_variants_after=\${n_var_after}
EOF

    echo "HWE filter: removed \${n_removed} variants (mode=\${filter_mode}, chrx=\${chrx_mode})"

    # ── HWE p-value distribution plot ─────────────────────────────────────────
    if command -v Rscript &>/dev/null; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
df <- read.table("hwe_stats.hwe", header=TRUE)
# Plot autosome distribution using the same test column used for filtering
test_col <- if (any(df$TEST == "UNAFF")) "UNAFF" else "ALL"
df_auto <- df[df$TEST == test_col & df$CHR < 23, ]
df_auto$neglog10p <- -log10(pmax(df_auto$P, 1e-300))
p <- ggplot(df_auto, aes(x=neglog10p)) +
    geom_histogram(bins=80, fill="steelblue", alpha=0.8) +
    geom_vline(xintercept=-log10(${params.hwe_p}), linetype="dashed", colour="red") +
    labs(title=paste0("HWE -log10(p): autosomes (test=", test_col, ")"),
         subtitle="Red line: autosome removal threshold",
         x="-log10(HWE p-value)", y="Number of variants") +
    theme_classic()
ggsave("hwe_plot.png", p, width=8, height=5)
RSCRIPT
    fi
    """
}
