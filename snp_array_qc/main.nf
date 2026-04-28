#!/usr/bin/env nextflow
/*
================================================================================
  SNP Array QC Pipeline — main.nf
================================================================================
  Workflow phases:
    01  Input validation
    02  Variant-level QC   (duplicate check, callrate, HWE, MAF)
    03  Sample-level QC    (missingness, sex, het, relatedness, PCA)
    04  Final report

  Variant-level QC can run per-chromosome for large datasets; sample-level QC
  always requires genome-wide data. Use --sample_qc_scope to control behaviour
  when only a chromosome subset is processed.

  Run:
    nextflow run snp_array_qc/main.nf \
      --bfile data/raw/genotypes \
      --outdir results/snp_array_qc \
      -profile singularity

  Chr-22 test run (variant QC only — sample QC set to provisional):
    nextflow run snp_array_qc/main.nf \
      --bfile data/raw/genotypes \
      --chroms 22 \
      --outdir results/snp_array_qc_chr22

  Skip a phase:
    nextflow run snp_array_qc/main.nf --bfile ... --run_sample_qc false
================================================================================
*/

nextflow.enable.dsl = 2

include { CHECK_VERSIONS          } from '../modules/check_versions'
include { INPUT_CHECK             } from './modules/input_check'
include { IMPUTATION_FILTER       } from './modules/imputation_filter'
include { DUPLICATE_VARIANT_CHECK } from './modules/duplicate_variant_check'
include { SAMPLE_CALLRATE         } from './modules/sample_callrate'
include { VARIANT_CALLRATE        } from './modules/variant_callrate'
include { SEX_CHECK               } from './modules/sex_check'
include { HETEROZYGOSITY          } from './modules/heterozygosity'
include { RELATEDNESS             } from './modules/relatedness'
include { ANCESTRY_PCA            } from './modules/ancestry_pca'
include { HWE_FILTER              } from './modules/hwe_filter'
include { MAF_FILTER              } from './modules/maf_filter'
include { MERGE_CHROMOSOMES       } from './modules/merge_chromosomes'
include { FINAL_REPORT            } from './modules/final_report'

// ── Chromosome range parser ───────────────────────────────────────────────────
// Accepts: "1-22", "22", "1,2,22", "all" (treated as 1-22)
def parseChroms(chrom_param) {
    def param = chrom_param.toString()
    if (param == "all") return (1..22).collect { it.toString() }
    if (param =~ /^\d+-\d+$/) {
        def (s, e) = param.split('-').collect { it.toInteger() }
        return (s..e).collect { it.toString() }
    }
    if (param =~ /,/) return param.split(',')*.trim()
    return [param.trim()]
}

// ── Determine effective sample QC scope ──────────────────────────────────────
def effectiveScope(chrom_param, declared_scope) {
    if (declared_scope != "auto") return declared_scope
    def chroms = parseChroms(chrom_param)
    def all22  = (1..22).every { n -> chroms.contains(n.toString()) || chroms.contains("chr${n}") }
    return all22 ? "genome_wide" : "provisional"
}

// ── Validate required parameters ─────────────────────────────────────────────
def validateParams() {
    if (!params.bfile) {
        error "ERROR: --bfile is required. Provide the PLINK binary prefix (without extension)."
    }
    if (!file("${params.bfile}.bed").exists()) error "ERROR: ${params.bfile}.bed not found"
    if (!file("${params.bfile}.bim").exists()) error "ERROR: ${params.bfile}.bim not found"
    if (!file("${params.bfile}.fam").exists()) error "ERROR: ${params.bfile}.fam not found"
    if (!["auto","genome_wide","provisional","skip"].contains(params.sample_qc_scope)) {
        error "ERROR: --sample_qc_scope must be one of: auto, genome_wide, provisional, skip"
    }
}

// ══════════════════════════════════════════════════════════════════════════════
//  MAIN WORKFLOW
// ══════════════════════════════════════════════════════════════════════════════
workflow {

    validateParams()

    def chrom_list = parseChroms(params.chroms)
    def scope      = effectiveScope(params.chroms, params.sample_qc_scope)

    CHECK_VERSIONS()

    log.info """
    ================================================================
    SNP Array QC Pipeline
    ================================================================
    PLINK prefix      : ${params.bfile}
    Output dir        : ${params.outdir}
    Phenotype file    : ${params.pheno ?: 'not provided'}
    Reference panel   : ${params.reference_panel ?: 'not provided'}
    HapMap info       : ${params.hapmap_info ?: 'not provided'}
    LD regions file   : ${params.ld_regions ?: 'not provided'}
    Chromosomes       : ${params.chroms}
    ----------------------------------------------------------------
    Phase switches
      Variant-level QC    : ${params.run_variant_qc}
      Sample-level QC     : ${params.run_sample_qc}
      Sample QC scope     : ${scope}${scope == 'provisional' ? '  *** PROVISIONAL — not all autosomes ***' : ''}
    ----------------------------------------------------------------
    Variant-level QC modules (Phase 2)
      Imputation filter   : ${params.run_imputation_filter}${params.run_imputation_filter ? "  (R2 >= ${params.imputation_r2}, file: ${params.info_file ?: 'none'})" : ''}
      Duplicate check     : ${params.run_duplicate_check}
      Variant missingness : ${params.run_variant_missingness}
      HWE filter          : ${params.run_hwe}
      MAF filter          : ${params.run_maf_filter}
    ----------------------------------------------------------------
    Sample-level QC modules (Phase 3)
      Sample missingness  : ${params.run_sample_missingness}
      Sex check           : ${params.run_sex_check}
      Heterozygosity      : ${params.run_heterozygosity}
      Relatedness         : ${params.run_relatedness}
      Ancestry PCA        : ${params.run_ancestry_pca}
    ----------------------------------------------------------------
    Thresholds
      Sample missingness  : ${params.sample_missingness}
      Variant missingness : ${params.variant_missingness}
      CC differential miss: ${params.cc_miss_p}
      HWE p (autosomes)   : ${params.hwe_p}
      HWE p (chrX)        : ${params.hwe_p_chrx}
      MAF                 : ${params.maf}
      Heterozygosity SD   : ${params.heterozygosity_sd}
      PI_HAT relatedness  : ${params.relatedness_pi_hat}
      PCA outlier SD      : ${params.pca_outlier_sd}
    ================================================================
    """.stripIndent()

    // ── Build input channels ──────────────────────────────────────────────────
    def prefix = file(params.bfile).name
    ch_plink = Channel.of([
        [ id: prefix ],
        file("${params.bfile}.bed"),
        file("${params.bfile}.bim"),
        file("${params.bfile}.fam")
    ])
    ch_pheno       = params.pheno           ? Channel.value(file(params.pheno))           : Channel.value([])
    ch_ref         = params.reference_panel ? Channel.value(params.reference_panel)        : Channel.value([])
    ch_ld_regions  = params.ld_regions      ? Channel.value(file(params.ld_regions))       : Channel.value([])
    ch_hapmap_info = params.hapmap_info     ? Channel.value(file(params.hapmap_info))      : Channel.value([])
    ch_info_file   = params.info_file       ? Channel.value(file(params.info_file))        : Channel.value([])

    ch_variant_exclusions = Channel.empty()
    ch_sample_exclusions  = Channel.empty()
    ch_qc_summaries       = Channel.empty()
    ch_qc_plots           = Channel.empty()
    ch_qc_data            = Channel.empty()

    // ════════════════════════════════════════════════════════════════════════
    //  PHASE 1 — Input validation
    // ════════════════════════════════════════════════════════════════════════
    INPUT_CHECK(ch_plink)
    ch_working      = INPUT_CHECK.out.plink
    ch_qc_summaries = ch_qc_summaries.mix(INPUT_CHECK.out.summary)

    // ════════════════════════════════════════════════════════════════════════
    //  PHASE 2 — Variant-level QC
    //
    //  Filters on variants/sites/genotypes. For SNP arrays the dataset is
    //  typically small enough to process genome-wide; the --chroms param
    //  and MERGE_CHROMOSOMES are provided for large multi-chip cohorts that
    //  benefit from chromosome-parallel execution.
    //
    //  Order: imputation filter → duplicates → callrate → HWE → MAF
    //  Imputation filter runs first to remove poorly imputed variants before
    //  any other QC metrics are computed. Disabled by default.
    //  HWE is placed before MAF so rare-variant HWE inflation doesn't bias
    //  the MAF cutoff; both run on the callrate-cleaned dataset.
    // ════════════════════════════════════════════════════════════════════════
    if (params.run_variant_qc) {

        if (params.run_imputation_filter) {
            IMPUTATION_FILTER(ch_working, ch_info_file)
            ch_working            = IMPUTATION_FILTER.out.plink
            ch_variant_exclusions = ch_variant_exclusions.mix(IMPUTATION_FILTER.out.removed_variants)
            ch_qc_summaries       = ch_qc_summaries.mix(IMPUTATION_FILTER.out.summary)
            ch_qc_plots           = ch_qc_plots.mix(IMPUTATION_FILTER.out.plots)
        } else {
            log.info "Imputation filter disabled (run_imputation_filter = false)"
        }

        if (params.run_duplicate_check) {
            DUPLICATE_VARIANT_CHECK(ch_working)
            ch_working            = DUPLICATE_VARIANT_CHECK.out.plink
            ch_variant_exclusions = ch_variant_exclusions.mix(DUPLICATE_VARIANT_CHECK.out.removed_variants)
            ch_qc_summaries       = ch_qc_summaries.mix(DUPLICATE_VARIANT_CHECK.out.summary)
        } else {
            log.warn "SKIPPING: Duplicate variant check (run_duplicate_check = false)"
        }

        if (params.run_variant_missingness) {
            VARIANT_CALLRATE(ch_working)
            ch_working            = VARIANT_CALLRATE.out.plink
            ch_variant_exclusions = ch_variant_exclusions.mix(VARIANT_CALLRATE.out.removed_variants)
            ch_qc_summaries       = ch_qc_summaries.mix(VARIANT_CALLRATE.out.summary)
            ch_qc_summaries       = ch_qc_summaries.mix(VARIANT_CALLRATE.out.cc_miss_summary)
            ch_qc_plots           = ch_qc_plots.mix(VARIANT_CALLRATE.out.plots)
        } else {
            log.warn "SKIPPING: Variant missingness (run_variant_missingness = false)"
        }

        if (params.run_hwe) {
            HWE_FILTER(ch_working, ch_pheno)
            ch_working            = HWE_FILTER.out.plink
            ch_variant_exclusions = ch_variant_exclusions.mix(HWE_FILTER.out.removed_variants)
            ch_qc_summaries       = ch_qc_summaries.mix(HWE_FILTER.out.summary)
            ch_qc_plots           = ch_qc_plots.mix(HWE_FILTER.out.plots)
        } else {
            log.warn "SKIPPING: HWE filter (run_hwe = false)"
        }

        if (params.run_maf_filter) {
            MAF_FILTER(ch_working)
            ch_working            = MAF_FILTER.out.plink
            ch_variant_exclusions = ch_variant_exclusions.mix(MAF_FILTER.out.removed_variants)
            ch_qc_summaries       = ch_qc_summaries.mix(MAF_FILTER.out.summary)
            ch_qc_plots           = ch_qc_plots.mix(MAF_FILTER.out.plots)
        } else {
            log.warn "SKIPPING: MAF filter (run_maf_filter = false)"
        }

        ch_all_excluded_variants = ch_variant_exclusions.collectFile(
            name:     'all_excluded_variants.txt',
            newLine:  true,
            storeDir: "${params.outdir}/exclusion_lists"
        )

    } else {
        log.warn "SKIPPING: Entire variant-level QC phase (run_variant_qc = false)"
        ch_all_excluded_variants = Channel.of('').collectFile(
            name: 'all_excluded_variants.txt',
            newLine: false,
            storeDir: "${params.outdir}/exclusion_lists"
        )
    }

    // ════════════════════════════════════════════════════════════════════════
    //  PHASE 3 — Sample-level QC
    //
    //  All modules here need genome-wide data to produce reliable metrics.
    //  When scope == "provisional" every module still runs, but the final
    //  report marks results as unsuitable for production filtering.
    //  When scope == "skip" the entire phase is omitted.
    // ════════════════════════════════════════════════════════════════════════
    if (params.run_sample_qc && scope != "skip") {

        if (scope == "provisional") {
            log.warn """
            *** PROVISIONAL SAMPLE QC ***
            Chromosomes processed : ${params.chroms}
            Not all autosomes were included — sample-level QC results are
            indicative only and must not be used for final filtering.
            Re-run with --chroms 1-22 for a genome-wide production run.
            """.stripIndent()
        }

        if (params.run_sample_missingness) {
            SAMPLE_CALLRATE(ch_working)
            ch_working           = SAMPLE_CALLRATE.out.plink
            ch_sample_exclusions = ch_sample_exclusions.mix(SAMPLE_CALLRATE.out.removed_samples)
            ch_qc_summaries      = ch_qc_summaries.mix(SAMPLE_CALLRATE.out.summary)
            ch_qc_data           = ch_qc_data.mix(SAMPLE_CALLRATE.out.imiss)
            ch_qc_data           = ch_qc_data.mix(SAMPLE_CALLRATE.out.removed_samples)
        } else {
            log.warn "SKIPPING: Sample missingness (run_sample_missingness = false)"
        }

        if (params.run_sex_check) {
            SEX_CHECK(ch_working, ch_pheno)
            ch_sample_exclusions = ch_sample_exclusions.mix(SEX_CHECK.out.discordant_samples)
            ch_qc_summaries      = ch_qc_summaries.mix(SEX_CHECK.out.summary)
            ch_qc_plots          = ch_qc_plots.mix(SEX_CHECK.out.plots)
            ch_qc_data           = ch_qc_data.mix(SEX_CHECK.out.sexcheck)
            ch_qc_data           = ch_qc_data.mix(SEX_CHECK.out.discordant_samples)
        } else {
            log.warn "SKIPPING: Sex check (run_sex_check = false)"
        }

        if (params.run_heterozygosity) {
            HETEROZYGOSITY(ch_working)
            ch_sample_exclusions = ch_sample_exclusions.mix(HETEROZYGOSITY.out.outlier_samples)
            ch_qc_summaries      = ch_qc_summaries.mix(HETEROZYGOSITY.out.summary)
            ch_qc_plots          = ch_qc_plots.mix(HETEROZYGOSITY.out.plots)
            ch_qc_data           = ch_qc_data.mix(HETEROZYGOSITY.out.het)
            ch_qc_data           = ch_qc_data.mix(HETEROZYGOSITY.out.outlier_samples)
        } else {
            log.warn "SKIPPING: Heterozygosity (run_heterozygosity = false)"
        }

        if (params.run_relatedness) {
            RELATEDNESS(ch_working, ch_ld_regions)
            ch_sample_exclusions = ch_sample_exclusions.mix(RELATEDNESS.out.related_samples)
            ch_qc_summaries      = ch_qc_summaries.mix(RELATEDNESS.out.summary)
            ch_qc_plots          = ch_qc_plots.mix(RELATEDNESS.out.plots)
            ch_qc_data           = ch_qc_data.mix(RELATEDNESS.out.related_samples)
        } else {
            log.warn "SKIPPING: Relatedness (run_relatedness = false)"
        }

        if (params.run_ancestry_pca) {
            ANCESTRY_PCA(ch_working, ch_ref, ch_ld_regions, ch_hapmap_info)
            ch_sample_exclusions = ch_sample_exclusions.mix(ANCESTRY_PCA.out.outlier_samples)
            ch_qc_summaries      = ch_qc_summaries.mix(ANCESTRY_PCA.out.summary)
            ch_qc_plots          = ch_qc_plots.mix(ANCESTRY_PCA.out.plots)
            ch_qc_data           = ch_qc_data.mix(ANCESTRY_PCA.out.outlier_samples)
        } else {
            log.warn "SKIPPING: Ancestry PCA (run_ancestry_pca = false)"
        }

        // Apply all sample exclusions in a single PLINK pass
        ch_all_excluded_samples = ch_sample_exclusions.collectFile(
            name:     'all_excluded_samples.txt',
            newLine:  true,
            storeDir: "${params.outdir}/exclusion_lists"
        )
        APPLY_SAMPLE_EXCLUSIONS(ch_working.combine(ch_all_excluded_samples))
        ch_final = APPLY_SAMPLE_EXCLUSIONS.out.plink

    } else {
        if (scope == "skip") {
            log.warn "SKIPPING: All sample-level QC (sample_qc_scope = skip)"
        } else {
            log.warn "SKIPPING: All sample-level QC (run_sample_qc = false)"
        }
        ch_final                = ch_working
        ch_all_excluded_samples = Channel.of('').collectFile(
            name: 'all_excluded_samples.txt',
            newLine: false,
            storeDir: "${params.outdir}/exclusion_lists"
        )
    }

    // ════════════════════════════════════════════════════════════════════════
    //  PHASE 4 — Final report
    // ════════════════════════════════════════════════════════════════════════
    if (params.run_final_report) {
        FINAL_REPORT(
            ch_final,
            ch_qc_summaries.collect(),
            ch_all_excluded_samples,
            ch_all_excluded_variants,
            ch_qc_plots.collect().ifEmpty([]),
            ch_qc_data.collect().ifEmpty([]),
            Channel.value(scope)
        )
    }
}

// ── Inline helper process: apply accumulated sample exclusions ────────────────
process APPLY_SAMPLE_EXCLUSIONS {
    label 'process_medium'
    publishDir "${params.outdir}/cleaned_data", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(exclusion_list)

    output:
    tuple val(meta), path("${meta.id}_sample_qc_pass.bed"),
                     path("${meta.id}_sample_qc_pass.bim"),
                     path("${meta.id}_sample_qc_pass.fam"), emit: plink

    script:
    """
    n_remove=\$(wc -l < ${exclusion_list} || echo 0)
    echo "Removing \${n_remove} samples flagged across all sample QC steps"

    plink \\
        --bfile ${bed.baseName} \\
        --remove ${exclusion_list} \\
        --make-bed \\
        --out ${meta.id}_sample_qc_pass \\
        --allow-no-sex

    echo "Samples retained: \$(wc -l < ${meta.id}_sample_qc_pass.fam)"
    """
}
