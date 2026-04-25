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

    script:
    def prefix = "${meta.id}"
    """
    n_var_before=\$(wc -l < ${bim})

    # ── Compute per-locus missingness, then filter ────────────────────────────
    plink \\
        --bfile ${bed.baseName} \\
        --missing \\
        --out variant_missingness \\
        --allow-no-sex

    plink \\
        --bfile ${bed.baseName} \\
        --geno ${params.variant_missingness} \\
        --make-bed \\
        --out ${prefix}_vmiss \\
        --allow-no-sex

    n_var_after=\$(wc -l < ${prefix}_vmiss.bim)
    n_removed=\$(( n_var_before - n_var_after ))

    # ── Extract IDs of removed variants ──────────────────────────────────────
    awk 'NR==FNR{seen[\$2]=1; next} !seen[\$2]{print \$2}' \\
        ${prefix}_vmiss.bim ${bim} > variant_callrate_removed.txt

    cat > variant_callrate_summary.txt << EOF
step=variant_callrate
dataset=${meta.id}
threshold=${params.variant_missingness}
n_variants_before=\${n_var_before}
n_variants_removed=\${n_removed}
n_variants_after=\${n_var_after}
EOF

    echo "Variant missingness filter: removed \${n_removed} of \${n_var_before} variants"
    """
}
