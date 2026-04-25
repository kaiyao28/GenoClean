/*
================================================================================
  MODULE: VARIANT_FILTERING
================================================================================
  Purpose:
    Apply hard filters or VQSR to a VCF to remove low-quality variants.
    Then apply per-genotype filters (GQ, DP) and recalculate missingness.

  Why this step:
    GATK Best Practices recommend filtering variants on quality annotation
    thresholds before analysis. Hard filtering is appropriate when:
      - Cohort is too small for VQSR (~<30 WES samples)
      - VQSR training resources are not available for the genome build
      - A rapid, transparent filter is preferred

    GATK recommends VQSR for larger cohorts with sufficient truth data.
    This module supports both modes via params.variant_filter_method.

  Hard-filter defaults (GATK recommendations):
    SNPs:   QD < 2.0, FS > 60, MQ < 40, MQRankSum < -12.5, ReadPosRS < -8
    indels: QD < 2.0, FS > 200, ReadPosRankSum < -20

  Genotype filters:
    GQ < params.min_gq  → set genotype to missing (./.)
    DP < params.min_dp  → set genotype to missing (./.)

  How to change any filter:
    nextflow run main.nf --snp_qd 3.0 --snp_fs 50.0

  How to disable:
    nextflow run main.nf --run_variant_filtering false

  Output:
    - filtered.vcf.gz         : hard-filtered VCF
    - genotype_filtered.vcf.gz : after per-genotype GQ/DP masking
    - variant_filtering_summary.txt : counts at each stage
================================================================================
*/

process VARIANT_FILTERING {
    label 'process_high'
    publishDir "${params.outdir}/cleaned_data",    mode: params.publish_dir_mode, pattern: "*.vcf.gz*"
    publishDir "${params.outdir}/qc_tables",       mode: params.publish_dir_mode, pattern: "*.txt"

    input:
    tuple val(meta), path(vcf)
    path reference_fasta

    output:
    tuple val(meta), path("${meta.id}_genotype_filtered.vcf.gz"),
                     path("${meta.id}_genotype_filtered.vcf.gz.tbi"), emit: vcf
    path "variant_filtering_summary.txt",                              emit: summary

    script:
    def prefix = "${meta.id}"
    def method = params.variant_filter_method
    """
    # ── Index input VCF if needed ──────────────────────────────────────────────
    if [ ! -f "${vcf}.tbi" ]; then
        bcftools index --tbi --threads ${task.cpus} ${vcf}
    fi

    n_before=\$(bcftools view --no-header ${vcf} | wc -l)

    ${method == "hard_filter" ? """
    # ── Hard filtering — SNPs ─────────────────────────────────────────────────
    # These thresholds follow GATK Best Practices for hard-filtering.
    # Each annotation is described in the GATK documentation:
    # QD  : QualByDepth — normalises QUAL by allele depth; low = unreliable
    # FS  : FisherStrand — strand bias; high = strand-specific artefact
    # MQ  : RMSMappingQuality — mapping quality of reads; low = alignment problem
    # MQRankSum : mapping quality difference between ref and alt reads
    # ReadPosRankSum : read position bias — variants near read ends are suspect
    gatk VariantFiltration \\
        -R ${reference_fasta} \\
        -V ${vcf} \\
        --filter-expression "QD < ${params.snp_qd}"               --filter-name "QD_lt_${params.snp_qd}" \\
        --filter-expression "FS > ${params.snp_fs}"               --filter-name "FS_gt_${params.snp_fs}" \\
        --filter-expression "MQ < ${params.snp_mq}"               --filter-name "MQ_lt_${params.snp_mq}" \\
        --filter-expression "MQRankSum < ${params.snp_mqrank_sum}" --filter-name "MQRankSum" \\
        --filter-expression "ReadPosRankSum < ${params.snp_readpos_rank_sum}" --filter-name "ReadPos" \\
        --filter-expression "QUAL < ${params.variant_qual}"        --filter-name "QUAL_lt_${params.variant_qual}" \\
        --genotype-filter-expression "GQ < ${params.min_gq}"      --genotype-filter-name "LowGQ" \\
        --genotype-filter-expression "DP < ${params.min_dp}"      --genotype-filter-name "LowDP" \\
        -O ${prefix}_site_filtered.vcf.gz

    # ── Additional indel filters ──────────────────────────────────────────────
    gatk VariantFiltration \\
        -R ${reference_fasta} \\
        -V ${prefix}_site_filtered.vcf.gz \\
        --filter-expression "VariantType == 'INDEL' && QD < ${params.indel_qd}"   --filter-name "INDEL_QD" \\
        --filter-expression "VariantType == 'INDEL' && FS > ${params.indel_fs}"   --filter-name "INDEL_FS" \\
        --filter-expression "VariantType == 'INDEL' && ReadPosRankSum < ${params.indel_readpos_rank_sum}" --filter-name "INDEL_ReadPos" \\
        -O ${prefix}_filtered.vcf.gz
    """ : """
    # ── VQSR filtering ────────────────────────────────────────────────────────
    # VQSR requires training resources (HapMap, OMNI, 1000G, dbSNP).
    # Ensure params.vqsr_hapmap, vqsr_omni, vqsr_1000g, vqsr_dbsnp are set.
    gatk VariantRecalibrator \\
        -V ${vcf} \\
        --trust-all-polymorphic \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15 ${params.vqsr_hapmap} \\
        --resource:omni,known=false,training=true,truth=true,prior=12 ${params.vqsr_omni} \\
        --resource:1000G,known=false,training=true,truth=false,prior=10 ${params.vqsr_1000g} \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${params.vqsr_dbsnp} \\
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \\
        -mode SNP \\
        -O snp_recal.vcf \\
        --tranches-file snp.tranches \\
        --rscript-file snp_plots.R

    gatk ApplyVQSR \\
        -V ${vcf} \\
        --recal-file snp_recal.vcf \\
        --tranches-file snp.tranches \\
        --truth-sensitivity-filter-level 99.5 \\
        --create-output-variant-index true \\
        -mode SNP \\
        -O ${prefix}_filtered.vcf.gz
    """}

    # ── Apply per-genotype filters (GQ and DP) ────────────────────────────────
    # Set low-quality genotypes to missing (./.) rather than removing the site.
    # This preserves variant sites that pass site-level filters.
    bcftools filter \\
        --threads ${task.cpus} \\
        --set-GTs . \\
        --include 'FMT/GQ >= ${params.min_gq} && FMT/DP >= ${params.min_dp}' \\
        ${prefix}_filtered.vcf.gz \\
        -O z -o ${prefix}_genotype_filtered.vcf.gz

    bcftools index --tbi --threads ${task.cpus} ${prefix}_genotype_filtered.vcf.gz

    # ── Count variants passing filters ────────────────────────────────────────
    n_pass=\$(bcftools view --no-header --apply-filters PASS ${prefix}_genotype_filtered.vcf.gz | wc -l)
    n_fail=\$(( n_before - n_pass ))

    cat > variant_filtering_summary.txt << EOF
step=variant_filtering
dataset=${meta.id}
filter_method=${method}
snp_qd_threshold=${params.snp_qd}
snp_fs_threshold=${params.snp_fs}
min_gq=${params.min_gq}
min_dp=${params.min_dp}
variant_qual=${params.variant_qual}
n_variants_before=\${n_before}
n_variants_filtered=\${n_fail}
n_variants_pass=\${n_pass}
EOF

    echo "Variant filtering: \${n_pass} variants PASS, \${n_fail} filtered"
    """
}
