/*
================================================================================
  MODULE: MERGE_CHROMOSOMES (SNP Array)
================================================================================
  Purpose:
    Merge chromosome-split PLINK binary datasets back into one genome-wide
    dataset before sample-level QC.
================================================================================
*/

process MERGE_CHROMOSOMES {
    label 'process_high'
    publishDir "${params.outdir}/cleaned_data", mode: params.publish_dir_mode

    input:
    path plink_files

    output:
    tuple val(["id": "merged"]), path("merged.bed"), path("merged.bim"), path("merged.fam"), emit: plink
    path "merge_chromosomes_summary.txt", emit: summary

    script:
    """
    ls *.bed | sed 's/.bed$//' | sort > merge_prefixes.txt
    head -n 1 merge_prefixes.txt > base_prefix.txt
    tail -n +2 merge_prefixes.txt > merge_list.txt
    base=\$(cat base_prefix.txt)

    if [ -s merge_list.txt ]; then
        plink --bfile "\${base}" --merge-list merge_list.txt --make-bed --out merged --allow-no-sex
    else
        cp "\${base}.bed" merged.bed
        cp "\${base}.bim" merged.bim
        cp "\${base}.fam" merged.fam
    fi

    cat > merge_chromosomes_summary.txt << EOF
step=merge_chromosomes
n_chromosome_files=\$(wc -l < merge_prefixes.txt)
n_samples=\$(wc -l < merged.fam)
n_variants=\$(wc -l < merged.bim)
EOF
    """
}
