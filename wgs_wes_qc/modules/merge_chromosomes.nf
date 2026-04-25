/*
================================================================================
  MODULE: MERGE_CHROMOSOMES (WGS/WES)
================================================================================
  Purpose:
    Concatenate chromosome-level VCFs back into a single genome-wide VCF before
    sample-level QC.
================================================================================
*/

process MERGE_CHROMOSOMES {
    label 'process_high'
    publishDir "${params.outdir}/cleaned_data", mode: params.publish_dir_mode

    input:
    path vcfs

    output:
    tuple val(["id": "merged", "input_type": params.input_type, "mode": params.mode]),
          path("merged.vcf.gz"),
          path("merged.vcf.gz.tbi"), emit: vcf
    path "merge_chromosomes_summary.txt", emit: summary

    script:
    """
    ls *.vcf.gz | sort > vcf_list.txt

    if [ \$(wc -l < vcf_list.txt) -gt 1 ]; then
        bcftools concat --threads ${task.cpus} -f vcf_list.txt -O z -o merged.vcf.gz
    else
        cp \$(cat vcf_list.txt) merged.vcf.gz
    fi

    bcftools index --tbi --threads ${task.cpus} merged.vcf.gz

    cat > merge_chromosomes_summary.txt << EOF
step=merge_chromosomes
n_chromosome_vcfs=\$(wc -l < vcf_list.txt)
n_variants=\$(bcftools view --no-header merged.vcf.gz | wc -l)
EOF
    """
}
