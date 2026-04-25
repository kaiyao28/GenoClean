/*
================================================================================
  MODULE: DUPLICATE_VARIANT_CHECK
================================================================================
  Purpose:
    Identify and remove duplicate SNP IDs from the PLINK dataset. Duplicates
    arise from multi-allelic sites split into separate entries, array design
    redundancy, or probe-name collisions. Retaining them causes problems in
    PLINK relatedness, PCA, and LD calculations.

  Why this step:
    PLINK will silently use the first occurrence of a duplicate ID, which can
    introduce subtle bias. Explicit removal is safer than relying on PLINK's
    handling of duplicates.

  Output:
    - Filtered PLINK files (deduplicated)
    - duplicate_variants.txt   : list of removed variant IDs
    - duplicate_summary.txt    : count of duplicates found

  How to disable:
    params.run_duplicate_check = false

  How to change behaviour:
    No threshold to change. If you want to KEEP the first occurrence instead
    of removing all copies, modify the awk logic below.
================================================================================
*/

process DUPLICATE_VARIANT_CHECK {
    label 'process_medium'
    publishDir "${params.outdir}/qc_tables", mode: params.publish_dir_mode, pattern: "*.txt"
    publishDir "${params.outdir}/cleaned_data", mode: params.publish_dir_mode, pattern: "*.{bed,bim,fam}", enabled: params.keep_intermediate

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("${meta.id}_dedup.bed"),
                     path("${meta.id}_dedup.bim"),
                     path("${meta.id}_dedup.fam"), emit: plink
    path "duplicate_variants.txt",                 emit: removed_variants
    path "duplicate_summary.txt",                  emit: summary

    script:
    def prefix = "${meta.id}"
    """
    # ── Find duplicate variant IDs (column 2 of .bim) ────────────────────────
    awk '{print \$2}' ${bim} | sort | uniq -d > duplicate_variants.txt

    n_dups=\$(wc -l < duplicate_variants.txt)

    cat > duplicate_summary.txt << EOF
step=duplicate_variant_check
dataset=${meta.id}
n_duplicates_removed=\${n_dups}
EOF

    if [ "\${n_dups}" -gt 0 ]; then
        echo "Removing \${n_dups} duplicate variant IDs"
        plink \\
            --bfile ${bed.baseName} \\
            --exclude duplicate_variants.txt \\
            --make-bed \\
            --out ${prefix}_dedup \\
            --allow-no-sex
    else
        echo "No duplicate variants found — copying files unchanged"
        # Create output files by running a no-op PLINK command
        plink \\
            --bfile ${bed.baseName} \\
            --make-bed \\
            --out ${prefix}_dedup \\
            --allow-no-sex
    fi
    """
}
