#!/usr/bin/env nextflow
/*
================================================================================
  WGS / WES QC Pipeline - main.nf
================================================================================
  Workflow phases:
    01  Input validation
    02  Variant-level QC   (chromosome-scoped VCF QC/filtering)
    03  Sample-level QC    (genome-wide BAM/CRAM and merged-VCF QC)
    04  Final report
================================================================================
*/

nextflow.enable.dsl = 2

include { CHECK_VERSIONS        } from '../modules/check_versions'
include { INPUT_CHECK           } from './modules/input_check'
include { FASTQC                } from './modules/fastqc'
include { ALIGNMENT_METRICS     } from './modules/alignment_metrics'
include { DUPLICATE_METRICS     } from './modules/duplicate_metrics'
include { COVERAGE_QC           } from './modules/coverage_qc'
include { CONTAMINATION_CHECK   } from './modules/contamination_check'
include { SEX_CHECK             } from './modules/sex_check'
include { VARIANT_CALLING_QC    } from './modules/variant_calling_qc'
include { VARIANT_FILTERING     } from './modules/variant_filtering'
include { SAMPLE_VARIANT_COUNTS } from './modules/sample_variant_counts'
include { RELATEDNESS           } from './modules/relatedness'
include { PCA_VARIANT_SELECTION } from './modules/pca_variant_selection'
include { ANCESTRY_PCA          } from './modules/ancestry_pca'
include { MERGE_CHROMOSOMES     } from './modules/merge_chromosomes'
include { FINAL_REPORT          } from './modules/final_report'

def parseChroms(chrom_param) {
    if (chrom_param == "all") return (1..22).collect { it.toString() }
    if (chrom_param =~ /^\d+-\d+$/) {
        def (s, e) = chrom_param.split('-').collect { it.toInteger() }
        return (s..e).collect { it.toString() }
    }
    if (chrom_param =~ /,/) return chrom_param.split(',')*.trim()
    return [chrom_param.trim()]
}

def effectiveScope(chrom_param, declared_scope) {
    if (declared_scope != "auto") return declared_scope
    def chroms = parseChroms(chrom_param)
    def all22 = (1..22).every { n -> chroms.contains(n.toString()) || chroms.contains("chr${n}") }
    return all22 ? "genome_wide" : "provisional"
}

def validateParams() {
    if (!params.samplesheet) error "ERROR: --samplesheet is required."
    if (!params.reference_fasta) error "ERROR: --reference_fasta is required."
    if (params.mode == "wes" && !params.target_intervals) {
        error "ERROR: --target_intervals is required when --mode wes"
    }
    if (!["fastq","bam","cram","vcf"].contains(params.input_type)) {
        error "ERROR: --input_type must be one of: fastq, bam, cram, vcf"
    }
    if (!["auto","genome_wide","provisional","skip"].contains(params.sample_qc_scope)) {
        error "ERROR: --sample_qc_scope must be one of: auto, genome_wide, provisional, skip"
    }
}

def parseSamplesheet(samplesheet_path) {
    def samples = []
    def lines = file(samplesheet_path).readLines()
    def header = lines[0].tokenize(',')*.trim()
    lines.drop(1).each { line ->
        def fields = line.tokenize(',')*.trim()
        def row = [header, fields].transpose().collectEntries()
        if (row.sample && (row.file1 || row.bam || row.fastq1 || row.vcf)) {
            samples << row
        }
    }
    return samples
}

workflow {
    validateParams()

    def chrom_list = parseChroms(params.chroms)
    def scope = effectiveScope(params.chroms, params.sample_qc_scope)

    CHECK_VERSIONS()

    log.info """
    ================================================================
    WGS / WES QC Pipeline
    ================================================================
    Mode              : ${params.mode.toUpperCase()}
    Input type        : ${params.input_type}
    Samplesheet       : ${params.samplesheet}
    Reference FASTA   : ${params.reference_fasta}
    Target intervals  : ${params.target_intervals ?: 'N/A'}
    Output dir        : ${params.outdir}
    Chromosomes       : ${params.chroms}
    ----------------------------------------------------------------
    Phase switches
      Variant-level QC    : ${params.run_variant_qc}
      Sample-level QC     : ${params.run_sample_qc}
      Sample QC scope     : ${scope}${scope == 'provisional' ? '  *** PROVISIONAL - not all autosomes ***' : ''}
    ----------------------------------------------------------------
    Variant-level modules
      Variant calling QC  : ${params.run_variant_calling_qc}
      Variant filtering   : ${params.run_variant_filtering}
    ----------------------------------------------------------------
    Sample-level modules
      FastQC              : ${params.run_fastqc}
      Alignment metrics   : ${params.run_alignment_metrics}
      Duplicate metrics   : ${params.run_duplicate_metrics}
      Coverage QC         : ${params.run_coverage_qc}
      Contamination       : ${params.run_contamination}
      Sex check           : ${params.run_sex_check_wgs}
      Sample counts       : ${params.run_sample_variant_counts}
      Relatedness         : ${params.run_relatedness_wgs}
      Ancestry PCA        : ${params.run_ancestry_pca_wgs}
    ================================================================
    """.stripIndent()

    ch_fasta = Channel.value(file(params.reference_fasta))
    ch_intervals = params.target_intervals ? Channel.value(file(params.target_intervals)) : Channel.value([])
    ch_ref_panel = params.reference_panel ? Channel.value(params.reference_panel) : Channel.value([])

    def samples = parseSamplesheet(params.samplesheet)
    ch_bam = Channel.empty()
    ch_fastq = Channel.empty()
    ch_vcf = Channel.empty()
    ch_input_check = Channel.empty()

    if (params.input_type in ["bam","cram"]) {
        ch_bam = Channel.fromList(samples).map { row ->
            def meta = [id: row.sample, input_type: params.input_type, mode: params.mode]
            def bam = file(row.file1 ?: row.bam)
            def bai = file("${bam}.bai").exists() ? file("${bam}.bai") : file("${bam.toString().replaceAll('\\.bam$','.bai')}")
            [meta, bam, bai]
        }
        ch_input_check = ch_bam
    } else if (params.input_type == "fastq") {
        ch_fastq = Channel.fromList(samples).map { row ->
            def meta = [id: row.sample, input_type: params.input_type, mode: params.mode]
            def r1 = file(row.file1 ?: row.fastq1)
            def r2 = (row.file2 ?: row.fastq2) ? file(row.file2 ?: row.fastq2) : []
            [meta, r1, r2]
        }
        ch_input_check = ch_fastq
    } else {
        ch_vcf = Channel.fromList(samples).map { row ->
            def meta = [id: row.sample, input_type: params.input_type, mode: params.mode]
            [meta, file(row.file1 ?: row.vcf)]
        }
        ch_input_check = ch_vcf.map { meta, vcf -> [meta, vcf, vcf] }
    }

    ch_qc_summaries = Channel.empty()
    ch_merged_vcf = Channel.empty()

    // PHASE 1 - Input validation
    INPUT_CHECK(ch_input_check, ch_fasta, ch_intervals)
    ch_qc_summaries = ch_qc_summaries.mix(INPUT_CHECK.out.summary)

    // PHASE 2 - Variant-level QC
    if (params.run_variant_qc && params.input_type == "vcf") {
        ch_vcf_by_chrom = ch_vcf.combine(Channel.fromList(chrom_list))
                            .map { meta, vcf, chrom -> [meta, vcf, chrom] }

        SELECT_CHROMOSOME(ch_vcf_by_chrom)
        ch_variant_working = SELECT_CHROMOSOME.out.vcf

        if (params.run_variant_calling_qc) {
            VARIANT_CALLING_QC(ch_variant_working, ch_fasta)
            ch_qc_summaries = ch_qc_summaries.mix(VARIANT_CALLING_QC.out.summary)
        } else {
            log.warn "SKIPPING: Variant calling QC (run_variant_calling_qc = false)"
        }

        if (params.run_variant_filtering) {
            VARIANT_FILTERING(ch_variant_working, ch_fasta)
            ch_variant_for_merge = VARIANT_FILTERING.out.vcf.map { meta, vcf, tbi -> vcf }
            ch_qc_summaries = ch_qc_summaries.mix(VARIANT_FILTERING.out.summary)
        } else {
            INDEX_CHROM_VCF(ch_variant_working)
            ch_variant_for_merge = INDEX_CHROM_VCF.out.vcf.map { meta, vcf, tbi -> vcf }
            ch_qc_summaries = ch_qc_summaries.mix(INDEX_CHROM_VCF.out.summary)
        }

        def merged_meta = [id: "merged", input_type: params.input_type, mode: params.mode]
        MERGE_CHROMOSOMES(ch_variant_for_merge.collect(), Channel.value(merged_meta))
        ch_merged_vcf = MERGE_CHROMOSOMES.out.vcf
        ch_qc_summaries = ch_qc_summaries.mix(MERGE_CHROMOSOMES.out.summary)
    } else {
        log.warn "SKIPPING: Variant-level QC phase (requires VCF input and run_variant_qc = true)"
        if (params.input_type == "vcf") {
            INDEX_INPUT_VCF(ch_vcf)
            ch_merged_vcf = INDEX_INPUT_VCF.out.vcf
            ch_qc_summaries = ch_qc_summaries.mix(INDEX_INPUT_VCF.out.summary)
        }
    }

    // PHASE 3 - Sample-level QC
    if (params.run_sample_qc && scope != "skip") {
        if (scope == "provisional") {
            log.warn """
            *** PROVISIONAL SAMPLE QC ***
            Chromosomes processed: ${params.chroms}
            Sample-level QC is not suitable for final filtering unless all
            autosomes are present. Re-run with --chroms 1-22 for production.
            """.stripIndent()
        }

        if (params.run_fastqc && params.input_type == "fastq") {
            FASTQC(ch_fastq)
            ch_qc_summaries = ch_qc_summaries.mix(FASTQC.out.summary)
        }

        if (params.run_alignment_metrics && params.input_type in ["bam","cram"]) {
            ALIGNMENT_METRICS(ch_bam, ch_fasta)
            ch_qc_summaries = ch_qc_summaries.mix(ALIGNMENT_METRICS.out.summary)
        }

        if (params.run_duplicate_metrics && params.input_type in ["bam","cram"]) {
            DUPLICATE_METRICS(ch_bam)
            ch_qc_summaries = ch_qc_summaries.mix(DUPLICATE_METRICS.out.summary)
        }

        if (params.run_coverage_qc && params.input_type in ["bam","cram"]) {
            COVERAGE_QC(ch_bam, ch_fasta, ch_intervals)
            ch_qc_summaries = ch_qc_summaries.mix(COVERAGE_QC.out.summary)
        }

        if (params.run_contamination && params.input_type in ["bam","cram"]) {
            CONTAMINATION_CHECK(ch_bam, ch_fasta)
            ch_qc_summaries = ch_qc_summaries.mix(CONTAMINATION_CHECK.out.summary)
        }

        if (params.run_sex_check_wgs && params.input_type in ["bam","cram"]) {
            SEX_CHECK(ch_bam, ch_fasta)
            ch_qc_summaries = ch_qc_summaries.mix(SEX_CHECK.out.summary)
        }

        if (params.input_type == "vcf") {
            if (params.run_sample_variant_counts) {
                SAMPLE_VARIANT_COUNTS(ch_merged_vcf)
                ch_qc_summaries = ch_qc_summaries.mix(SAMPLE_VARIANT_COUNTS.out.summary)
            }

            if (params.run_relatedness_wgs) {
                RELATEDNESS(ch_merged_vcf)
                ch_qc_summaries = ch_qc_summaries.mix(RELATEDNESS.out.summary)
            }

            if (params.run_ancestry_pca_wgs) {
                PCA_VARIANT_SELECTION(ch_merged_vcf)
                ch_qc_summaries = ch_qc_summaries.mix(PCA_VARIANT_SELECTION.out.summary)

                ANCESTRY_PCA(PCA_VARIANT_SELECTION.out.vcf, ch_ref_panel)
                ch_qc_summaries = ch_qc_summaries.mix(ANCESTRY_PCA.out.summary)
            }
        }
    } else {
        log.warn "SKIPPING: Sample-level QC phase"
    }

    // PHASE 4 - Final report
    if (params.run_final_report) {
        FINAL_REPORT(ch_qc_summaries.collect())
    }
}

process SELECT_CHROMOSOME {
    label 'process_medium'
    publishDir "${params.outdir}/variant_qc/chromosomes", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), val(chrom)

    output:
    tuple val(meta + [id: "${meta.id}.chr${chrom}", chrom: chrom]), path("${meta.id}.chr${chrom}.vcf.gz"), emit: vcf

    script:
    """
    if [[ "${vcf}" == *.vcf.gz ]]; then
        cp ${vcf} input.vcf.gz
    else
        bgzip -c ${vcf} > input.vcf.gz
    fi
    bcftools index --tbi --force --threads ${task.cpus} input.vcf.gz

    region1="${chrom}"
    region2="chr${chrom}"

    if bcftools view --regions "\${region1}" --no-header input.vcf.gz | head -n 1 | grep -q .; then
        region="\${region1}"
    else
        region="\${region2}"
    fi

    bcftools view --threads ${task.cpus} --regions "\${region}" input.vcf.gz -O z -o ${meta.id}.chr${chrom}.vcf.gz
    bcftools index --tbi --threads ${task.cpus} ${meta.id}.chr${chrom}.vcf.gz
    """
}

process INDEX_CHROM_VCF {
    label 'process_low'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.indexed.vcf.gz"), path("${meta.id}.indexed.vcf.gz.tbi"), emit: vcf
    path "index_chrom_vcf_summary.txt", emit: summary

    script:
    """
    bcftools view --threads ${task.cpus} ${vcf} -O z -o ${meta.id}.indexed.vcf.gz
    bcftools index --tbi --threads ${task.cpus} ${meta.id}.indexed.vcf.gz

    cat > index_chrom_vcf_summary.txt << EOF
step=index_chrom_vcf
dataset=${meta.id}
n_variants=\$(bcftools view --no-header ${meta.id}.indexed.vcf.gz | wc -l)
EOF
    """
}

process INDEX_INPUT_VCF {
    label 'process_low'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.indexed.vcf.gz"), path("${meta.id}.indexed.vcf.gz.tbi"), emit: vcf
    path "index_input_vcf_summary.txt", emit: summary

    script:
    """
    bcftools view --threads ${task.cpus} ${vcf} -O z -o ${meta.id}.indexed.vcf.gz
    bcftools index --tbi --threads ${task.cpus} ${meta.id}.indexed.vcf.gz

    cat > index_input_vcf_summary.txt << EOF
step=index_input_vcf
dataset=${meta.id}
n_variants=\$(bcftools view --no-header ${meta.id}.indexed.vcf.gz | wc -l)
EOF
    """
}
