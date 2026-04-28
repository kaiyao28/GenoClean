/*
================================================================================
  MODULE: VARIANT_CALLRATE
================================================================================
  Purpose:
    Remove variants (SNPs) with a high proportion of missing calls across
    samples. Poorly genotyped variants are unreliable for association testing
    and inflate test statistics.

  Why this step:
    A two-pass strategy is common: a loose filter (e.g. 10%) before sample QC
    to avoid sample exclusion being driven by a few bad SNPs, and a tighter
    filter (e.g. 2%) after sample QC. This module implements the tighter pass.
    Anderson et al. 2010 recommend 2% as a conservative standard threshold.

  Default threshold: params.variant_missingness = 0.02 (2%)

  How to change:
    nextflow run main.nf --variant_missingness 0.05

  How to disable:
    nextflow run main.nf --run_variant_missingness false

  Output:
    - Filtered PLINK files (variants passing QC)
    - variant_callrate_removed.txt   : SNP IDs of removed variants
    - variant_callrate_summary.txt   : counts before/after
    - variant_missingness.lmiss      : per-variant missingness table
================================================================================
*/

process VARIANT_CALLRATE {
    label 'process_medium'
    publishDir "${params.outdir}/qc_tables",   mode: params.publish_dir_mode, pattern: "*.{lmiss,txt}"
    publishDir "${params.outdir}/qc_plots",   mode: params.publish_dir_mode, pattern: "*.png"
    publishDir "${params.outdir}/cleaned_data", mode: params.publish_dir_mode, pattern: "*.{bed,bim,fam}", enabled: params.keep_intermediate

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("${meta.id}_vmiss.bed"),
                     path("${meta.id}_vmiss.bim"),
                     path("${meta.id}_vmiss.fam"), emit: plink
    path "variant_callrate_removed.txt",            emit: removed_variants
    path "variant_callrate_summary.txt",            emit: summary
    path "variant_missingness.lmiss",               emit: lmiss
    path "cc_miss_removed.txt",                     optional: true, emit: cc_miss_removed
    path "cc_miss_summary.txt",                     optional: true, emit: cc_miss_summary
    path "*.png",                                   optional: true, emit: plots

    script:
    def prefix = "${meta.id}"
    """
    n_var_before=\$(wc -l < ${bim})

    # ── Step 1: Compute per-locus missingness stats ───────────────────────────
    plink \\
        --bfile ${bed.baseName} \\
        --missing \\
        --out variant_missingness \\
        --allow-no-sex

    # ── Step 2: Case-control differential missingness test ────────────────────
    # Run on pre-filter data so batch-specific call-rate differences in variants
    # that pass the overall threshold are not missed. Anderson et al. 2010.
    n_cases=\$(awk '\$6==2' ${fam} | wc -l)
    n_controls=\$(awk '\$6==1' ${fam} | wc -l)

    touch cc_miss_removed.txt
    touch cc_miss_summary.txt

    if [ "\${n_cases}" -gt 0 ] && [ "\${n_controls}" -gt 0 ]; then
        echo "Case-control data detected — running differential missingness test"
        plink \\
            --bfile ${bed.baseName} \\
            --test-missing \\
            --out cc_miss \\
            --allow-no-sex

        awk -v p="${params.cc_miss_p}" 'NR>1 && \$5+0 < p {print \$2}' \\
            cc_miss.missing > cc_miss_removed.txt

        n_cc_fail=\$(wc -l < cc_miss_removed.txt)

        cat > cc_miss_summary.txt << EOF
step=cc_miss_test
dataset=${meta.id}
cc_miss_p_threshold=${params.cc_miss_p}
n_cases=\${n_cases}
n_controls=\${n_controls}
n_variants_failing=\${n_cc_fail}
EOF
        echo "Case-control missingness: \${n_cc_fail} variants flagged (p < ${params.cc_miss_p})"

        if command -v Rscript &>/dev/null && [ -f cc_miss.missing ]; then
            Rscript - << 'RSCRIPT'
library(ggplot2)
df <- read.table("cc_miss.missing", header=TRUE)
df$neglog10p <- -log10(pmax(df$P, 1e-300))
p <- ggplot(df, aes(x=neglog10p)) +
    geom_histogram(bins=80, fill="steelblue", alpha=0.8) +
    geom_vline(xintercept=-log10(${params.cc_miss_p}), linetype="dashed", colour="red") +
    labs(title="Case-control differential missingness: -log10(p)",
         subtitle="Red line: removal threshold — variants to the right are flagged",
         x="-log10(p-value)", y="Number of variants") +
    theme_classic()
ggsave("cc_miss_plot.png", p, width=8, height=5)
RSCRIPT
        fi
    else
        echo "No case-control phenotype coding — skipping differential missingness test"
    fi

    # ── Step 3: Apply overall missingness filter + cc_miss exclusion ──────────
    # Both filters are applied in a single PLINK call so the output files reflect
    # all variant-level removals performed in this step.
    plink \\
        --bfile ${bed.baseName} \\
        --geno ${params.variant_missingness} \\
        --exclude cc_miss_removed.txt \\
        --make-bed \\
        --out ${prefix}_vmiss \\
        --allow-no-sex

    n_var_after=\$(wc -l < ${prefix}_vmiss.bim)
    n_removed=\$(( n_var_before - n_var_after ))

    # ── Extract IDs of all removed variants ──────────────────────────────────
    awk 'NR==FNR{seen[\$2]=1; next} !seen[\$2]{print \$2}' \\
        ${prefix}_vmiss.bim ${bim} > variant_callrate_removed.txt

    cat > variant_callrate_summary.txt << EOF
step=variant_callrate
dataset=${meta.id}
threshold=${params.variant_missingness}
cc_miss_p=${params.cc_miss_p}
n_variants_before=\${n_var_before}
n_variants_removed=\${n_removed}
n_variants_after=\${n_var_after}
EOF

    echo "Variant missingness filter: removed \${n_removed} of \${n_var_before} variants"

    # ── Variant missingness distribution plot ─────────────────────────────────
    if command -v Rscript &>/dev/null; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
df <- read.table("variant_missingness.lmiss", header=TRUE)
p <- ggplot(df, aes(x=F_MISS)) +
    geom_histogram(bins=80, fill="steelblue", alpha=0.8) +
    geom_vline(xintercept=${params.variant_missingness}, linetype="dashed", colour="red") +
    labs(title="Variant missingness distribution",
         subtitle="Red line: removal threshold",
         x="Proportion of missing calls per variant", y="Number of variants") +
    theme_classic()
ggsave("vmiss_plot.png", p, width=8, height=5)
RSCRIPT
    fi
    """
}
