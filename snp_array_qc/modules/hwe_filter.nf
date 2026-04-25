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

    script:
    def prefix = "${meta.id}"
    """
    n_var_before=\$(wc -l < ${bim})

    # ── Detect whether case/control phenotypes are coded in the FAM file ──────
    # FAM column 6: 1=control, 2=case; 0 or -9 = missing
    n_cases=\$(awk '\$6==2' ${fam} | wc -l)
    n_controls=\$(awk '\$6==1' ${fam} | wc -l)

    if [ "\${n_controls}" -gt 0 ]; then
        echo "Case-control data detected (\${n_cases} cases, \${n_controls} controls)"
        echo "Applying HWE filter in CONTROLS ONLY (as recommended)"
        # --hwe filters on controls when phenotype is coded; --hwe-all overrides this
        plink \\
            --bfile ${bed.baseName} \\
            --hwe ${params.hwe_p} \\
            --make-bed \\
            --out ${prefix}_hwe \\
            --allow-no-sex
        filter_mode="controls_only"
    else
        echo "No case-control coding detected — applying HWE filter across all samples"
        plink \\
            --bfile ${bed.baseName} \\
            --hwe ${params.hwe_p} \\
            --make-bed \\
            --out ${prefix}_hwe \\
            --allow-no-sex
        filter_mode="all_samples"
    fi

    n_var_after=\$(wc -l < ${prefix}_hwe.bim)
    n_removed=\$(( n_var_before - n_var_after ))

    # ── Record removed variant IDs ────────────────────────────────────────────
    awk 'NR==FNR{seen[\$2]=1; next} !seen[\$2]{print \$2}' \\
        ${prefix}_hwe.bim ${bim} > hwe_removed_variants.txt

    cat > hwe_summary.txt << EOF
step=hwe_filter
dataset=${meta.id}
hwe_p_threshold=${params.hwe_p}
filter_mode=\${filter_mode}
n_variants_before=\${n_var_before}
n_variants_removed=\${n_removed}
n_variants_after=\${n_var_after}
EOF

    echo "HWE filter: removed \${n_removed} variants (p < ${params.hwe_p})"
    """
}
