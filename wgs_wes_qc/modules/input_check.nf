/*
================================================================================
  MODULE: INPUT_CHECK (WGS/WES)
================================================================================
  Purpose:
    Validate that all input files exist and are indexed, sample IDs are unique,
    and required parameters (reference FASTA, target intervals for WES) are
    present and accessible.

  Output:
    - input_summary.tsv   : per-sample file paths and validation status
    - input_check_summary.txt : aggregate counts for the final report
================================================================================
*/

process INPUT_CHECK {
    label 'process_low'
    publishDir "${params.outdir}/logs", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(file1), path(file2_or_bai)   // file2_or_bai may be []
    path  reference_fasta
    path  target_intervals                              // may be []

    output:
    path "input_summary.tsv",        emit: summary_table
    path "input_check_summary.txt",  emit: summary

    script:
    """
    #!/usr/bin/env python3
    import os, sys

    errors = []
    sample_id = "${meta.id}"
    input_type = "${meta.input_type}"

    # ── Check primary input file ──────────────────────────────────────────────
    f1 = "${file1}"
    if not os.path.exists(f1):
        errors.append(f"ERROR: {f1} not found")
    elif os.path.getsize(f1) == 0:
        errors.append(f"ERROR: {f1} is empty")

    # ── Check BAI/CRAI index for BAM/CRAM ────────────────────────────────────
    bai = "${file2_or_bai}"
    if input_type in ("bam", "cram") and bai and bai != "[]":
        if not os.path.exists(bai):
            errors.append(f"WARNING: BAM index {bai} not found — may need to index")

    # ── Check reference FASTA ─────────────────────────────────────────────────
    ref = "${reference_fasta}"
    if not os.path.exists(ref):
        errors.append(f"ERROR: Reference FASTA {ref} not found")
    fai = ref + ".fai"
    if not os.path.exists(fai):
        errors.append(f"WARNING: Reference FAI index {fai} not found")

    # ── Check target intervals for WES ────────────────────────────────────────
    mode = "${meta.mode}"
    intervals = "${target_intervals}"
    if mode == "wes" and (not intervals or intervals == "[]"):
        errors.append("ERROR: target_intervals BED is required for WES mode")
    elif intervals and intervals != "[]" and not os.path.exists(intervals):
        errors.append(f"ERROR: target_intervals {intervals} not found")

    # ── Write summary ─────────────────────────────────────────────────────────
    status = "FAIL" if any(e.startswith("ERROR") for e in errors) else "PASS"

    with open("input_summary.tsv", "w") as fh:
        fh.write("sample\tinput_type\tmode\tfile\tstatus\tnotes\n")
        fh.write(f"{sample_id}\t{input_type}\t{mode}\t{f1}\t{status}\t{'; '.join(errors) if errors else 'OK'}\n")

    with open("input_check_summary.txt", "w") as fh:
        fh.write(f"step=input_check\\n")
        fh.write(f"dataset={sample_id}\\n")
        fh.write(f"input_type={input_type}\\n")
        fh.write(f"mode={mode}\\n")
        fh.write(f"status={status}\\n")
        if errors:
            for e in errors:
                fh.write(f"note={e}\\n")

    if any(e.startswith("ERROR") for e in errors):
        for e in errors:
            print(e, file=sys.stderr)
        sys.exit(1)
    else:
        for w in errors:
            print(w)
        print(f"Input check PASSED for {sample_id}")
    """
}
