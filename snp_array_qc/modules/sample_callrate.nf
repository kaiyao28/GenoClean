/*
================================================================================
  MODULE: SAMPLE_CALLRATE
================================================================================
  Purpose:
    Remove individuals with a high proportion of missing genotype calls.
    High missingness in a sample typically reflects poor DNA quality,
    hybridisation failure, or sample degradation.

  Why this step:
    Anderson et al. 2010 (Nature Protocols) recommend filtering samples with
    > 2% missing genotypes as a first-pass quality gate. Samples with many
    missing genotypes also distort downstream statistics such as relatedness
    estimation and PCA.

  Default threshold: params.sample_missingness = 0.02 (2%)

  How to change:
    nextflow run main.nf --sample_missingness 0.05

  How to disable:
    nextflow run main.nf --run_sample_missingness false

  Output:
    - Filtered PLINK files (samples passing QC)
    - sample_callrate_removed.txt   : FID IID of removed samples
    - sample_callrate_summary.txt   : counts before/after
    - sample_missingness.imiss      : per-sample missingness table
================================================================================
*/

process SAMPLE_CALLRATE {
    label 'process_medium'
    publishDir "${params.outdir}/qc_tables",   mode: params.publish_dir_mode, pattern: "*.{imiss,txt}"
    publishDir "${params.outdir}/cleaned_data", mode: params.publish_dir_mode, pattern: "*.{bed,bim,fam}", enabled: params.keep_intermediate

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("${meta.id}_smiss.bed"),
                     path("${meta.id}_smiss.bim"),
                     path("${meta.id}_smiss.fam"), emit: plink
    path "sample_callrate_removed.txt",             emit: removed_samples
    path "sample_callrate_summary.txt",             emit: summary
    path "sample_missingness.imiss",                emit: imiss

    script:
    def prefix = "${meta.id}"
    """
    n_before=\$(wc -l < ${fam})

    # ── Compute per-sample missingness, then filter ───────────────────────────
    # --mind removes samples with missingness > threshold
    # --missing writes the .imiss file for downstream plotting
    plink \\
        --bfile ${bed.baseName} \\
        --missing \\
        --out sample_missingness \\
        --allow-no-sex

    plink \\
        --bfile ${bed.baseName} \\
        --mind ${params.sample_missingness} \\
        --make-bed \\
        --out ${prefix}_smiss \\
        --allow-no-sex

    n_after=\$(wc -l < ${prefix}_smiss.fam)
    n_removed=\$(( n_before - n_after ))

    # ── Extract FID/IID of removed samples ───────────────────────────────────
    # Compare FAM files to get removed individuals
    awk 'NR==FNR{seen[\$1"\t"\$2]=1; next} !seen[\$1"\t"\$2]{print \$1"\t"\$2}' \\
        ${prefix}_smiss.fam ${fam} > sample_callrate_removed.txt

    cat > sample_callrate_summary.txt << EOF
step=sample_callrate
dataset=${meta.id}
threshold=${params.sample_missingness}
n_before=\${n_before}
n_removed=\${n_removed}
n_after=\${n_after}
EOF

    echo "Sample missingness filter: removed \${n_removed} of \${n_before} samples"
    """
}
