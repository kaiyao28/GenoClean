/*
================================================================================
  MODULE: INPUT_CHECK (WGS/WES)
================================================================================
  Purpose:
    Validate that input files and required reference files are present before
    downstream QC modules run.

  Output:
    - input_summary.tsv
    - input_check_summary.txt
================================================================================
*/

process INPUT_CHECK {
    label 'process_low'
    publishDir "${params.outdir}/logs", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(file1), path(file2_or_bai)
    path reference_fasta
    path target_intervals

    output:
    path "input_summary.tsv",       emit: summary_table
    path "input_check_summary.txt", emit: summary

    script:
    """
    notes=""
    status="PASS"

    if [ ! -s "${file1}" ]; then
        notes="\${notes}ERROR: ${file1} missing or empty; "
        status="FAIL"
    fi

    if [ "${meta.input_type}" = "bam" ] || [ "${meta.input_type}" = "cram" ]; then
        if [ "${file2_or_bai}" != "[]" ] && [ ! -e "${file2_or_bai}" ]; then
            notes="\${notes}WARNING: index ${file2_or_bai} not found; "
        fi
    fi

    if [ ! -s "${reference_fasta}" ]; then
        notes="\${notes}ERROR: reference FASTA ${reference_fasta} missing or empty; "
        status="FAIL"
    fi

    if [ ! -e "${reference_fasta}.fai" ]; then
        notes="\${notes}WARNING: reference FASTA index ${reference_fasta}.fai not found; "
    fi

    if [ "${meta.mode}" = "wes" ]; then
        if [ "${target_intervals}" = "[]" ] || [ ! -s "${target_intervals}" ]; then
            notes="\${notes}ERROR: target intervals required for WES mode; "
            status="FAIL"
        fi
    elif [ -n "${target_intervals}" ] && [ "${target_intervals}" != "[]" ] && [ ! -s "${target_intervals}" ]; then
        notes="\${notes}ERROR: target intervals ${target_intervals} not found; "
        status="FAIL"
    fi

    if [ -z "\${notes}" ]; then
        notes="OK"
    fi

    cat > input_summary.tsv << EOF
sample	input_type	mode	file	status	notes
${meta.id}	${meta.input_type}	${meta.mode}	${file1}	\${status}	\${notes}
EOF

    cat > input_check_summary.txt << EOF
step=input_check
dataset=${meta.id}
input_type=${meta.input_type}
mode=${meta.mode}
status=\${status}
note=\${notes}
EOF

    if [ "\${status}" = "FAIL" ]; then
        echo "\${notes}" >&2
        exit 1
    fi

    echo "Input check \${status} for ${meta.id}: \${notes}"
    """
}
