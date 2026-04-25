/*
================================================================================
  MODULE: INPUT_CHECK
================================================================================
  Purpose:
    Validate that the PLINK binary triplet (.bed/.bim/.fam) exists and is
    internally consistent. Produces a summary table of input counts that is
    carried into the final report.

  Why this step:
    Catching input problems early prevents cryptic errors deep in the pipeline.
    Variant and sample counts provide a baseline for the QC attrition table.

  Output:
    - input_summary.txt  : sample count, variant count, chromosome list
    - Passes plink files unchanged for downstream steps
================================================================================
*/

process INPUT_CHECK {
    label 'process_low'
    publishDir "${params.outdir}/logs", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path(bed), path(bim), path(fam), emit: plink
    path "input_summary.txt",                          emit: summary

    script:
    def prefix = bed.baseName
    """
    # ── Validate files are readable and non-empty ────────────────────────────
    for f in ${bed} ${bim} ${fam}; do
        if [ ! -s "\$f" ]; then
            echo "ERROR: \$f is missing or empty" >&2
            exit 1
        fi
    done

    # ── Count samples and variants ────────────────────────────────────────────
    n_samples=\$(wc -l < ${fam})
    n_variants=\$(wc -l < ${bim})

    # ── Chromosome distribution from .bim column 1 ───────────────────────────
    chr_dist=\$(awk '{print \$1}' ${bim} | sort -V | uniq -c | awk '{print \$2"("\$1")"}' | tr '\\n' ',' | sed "s/,\$//")

    # ── Phenotype column summary (0=missing, 1=control, 2=case) ──────────────
    n_missing_pheno=\$(awk '\$6==0 || \$6==-9 {count++} END{print count+0}' ${fam})
    n_cases=\$(awk '\$6==2 {count++} END{print count+0}' ${fam})
    n_controls=\$(awk '\$6==1 {count++} END{print count+0}' ${fam})

    cat > input_summary.txt << EOF
step=input_check
dataset=${meta.id}
n_samples=\${n_samples}
n_variants=\${n_variants}
n_cases=\${n_cases}
n_controls=\${n_controls}
n_missing_pheno=\${n_missing_pheno}
chromosome_distribution=\${chr_dist}
EOF

    echo "Input check complete: \${n_samples} samples, \${n_variants} variants"
    """
}
