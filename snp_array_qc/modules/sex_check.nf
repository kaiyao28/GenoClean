/*
================================================================================
  MODULE: SEX_CHECK
================================================================================
  Purpose:
    Compare the genetically inferred sex (based on X-chromosome homozygosity)
    with the reported sex in the FAM file. Discordant samples likely represent
    sample swaps, mislabelling, or pseudoautosomal region artefacts.

  Why this step:
    Sex discordance is a strong indicator of sample-level error. These samples
    should be excluded before any sex-stratified analysis. PLINK calculates the
    X-chromosome inbreeding coefficient (F statistic): females (two X) have
    F ≈ 0; males (one X) have F ≈ 1.

  Thresholds:
    params.sex_check_f_lower_female = 0.2   (F < 0.2 → inferred female)
    params.sex_check_f_upper_male   = 0.8   (F > 0.8 → inferred male)
    Samples with 0.2 ≤ F ≤ 0.8 have ambiguous inferred sex and are flagged.

  How to change:
    nextflow run main.nf --sex_check_f_lower_female 0.3 --sex_check_f_upper_male 0.7

  How to disable:
    nextflow run main.nf --run_sex_check false

  Output:
    - sex_check.sexcheck      : PLINK sex-check output (all samples)
    - sex_discordant.txt      : FID IID of samples with sex mismatch
    - sex_check_summary.txt   : counts of concordant / discordant / ambiguous
================================================================================
*/

process SEX_CHECK {
    label 'process_medium'
    publishDir "${params.outdir}/qc_tables", mode: params.publish_dir_mode, pattern: "*.{sexcheck,txt}"
    publishDir "${params.outdir}/qc_plots",  mode: params.publish_dir_mode, pattern: "*.png"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path pheno   // may be empty list if not provided

    output:
    path "sex_discordant.txt",      emit: discordant_samples
    path "sex_check_summary.txt",   emit: summary
    path "sex_check.sexcheck",      emit: sexcheck
    path "*.png",                   optional: true, emit: plots

    script:
    def prefix = "${meta.id}"
    """
    # ── Run PLINK sex check ───────────────────────────────────────────────────
    # PLINK uses X-chromosome SNPs to estimate the F inbreeding coefficient.
    # Requires chrX SNPs to be present; the process will warn if absent.
    plink \\
        --bfile ${bed.baseName} \\
        --check-sex ${params.sex_check_f_lower_female} ${params.sex_check_f_upper_male} \\
        --out sex_check \\
        --allow-no-sex

    # ── Parse sex_check.sexcheck for discordant and ambiguous samples ─────────
    # Columns: FID IID PEDSEX SNPSEX STATUS F
    # PEDSEX: 1=male 2=female 0=unknown; SNPSEX: same encoding
    # STATUS: OK or PROBLEM (discordant or ambiguous)

    awk 'NR>1 && \$5=="PROBLEM" {print \$1"\t"\$2}' sex_check.sexcheck > sex_discordant.txt

    n_total=\$(awk 'NR>1' sex_check.sexcheck | wc -l)
    n_problem=\$(wc -l < sex_discordant.txt)
    n_ok=\$(( n_total - n_problem ))

    n_ambiguous=\$(awk 'NR>1 && \$4==0 && \$5=="PROBLEM"' sex_check.sexcheck | wc -l)
    n_discordant=\$(( n_problem - n_ambiguous ))

    cat > sex_check_summary.txt << EOF
step=sex_check
dataset=${meta.id}
f_lower_female=${params.sex_check_f_lower_female}
f_upper_male=${params.sex_check_f_upper_male}
n_total=\${n_total}
n_concordant=\${n_ok}
n_discordant=\${n_discordant}
n_ambiguous=\${n_ambiguous}
n_flagged=\${n_problem}
EOF

    echo "Sex check: \${n_problem} samples flagged (discordant or ambiguous)"

    # ── Optional: generate F-statistic plot if Rscript is available ──────────
    if command -v Rscript &>/dev/null; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
df <- read.table("sex_check.sexcheck", header=TRUE)
df\$sex_label <- ifelse(df\$PEDSEX==1, "Reported male",
                  ifelse(df\$PEDSEX==2, "Reported female", "Unknown"))
p <- ggplot(df, aes(x=F, fill=sex_label)) +
    geom_histogram(bins=80, alpha=0.7) +
    geom_vline(xintercept=c(${params.sex_check_f_lower_female}, ${params.sex_check_f_upper_male}),
               linetype="dashed", colour="red") +
    labs(title="Sex check: X-chromosome F statistic",
         x="F statistic (0=female, 1=male)", y="Count",
         fill="Reported sex") +
    theme_classic()
ggsave("sex_check_F_stat.png", p, width=8, height=5)
RSCRIPT
    fi
    """
}
