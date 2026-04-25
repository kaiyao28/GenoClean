/*
================================================================================
  MODULE: DUPLICATE_METRICS
================================================================================
  Purpose:
    Mark PCR/optical duplicates using Picard MarkDuplicates and report the
    duplication rate. Samples with very high duplication rates may have
    insufficient unique molecules for reliable variant calling.

  Why this step:
    PCR duplicates inflate the apparent depth without adding information and
    can introduce systematic errors in variant allele frequencies. For WES,
    high duplication is common with tight capture kits; for WGS it typically
    indicates low input DNA or over-amplification.

  Default threshold: params.max_duplication_rate = 0.20 (20%)
  Samples exceeding this threshold are flagged but NOT automatically removed,
  because the threshold may need to be adjusted per-experiment.

  How to change:
    params.max_duplication_rate = 0.30  # allow 30% for WES

  How to disable:
    params.run_duplicate_metrics = false

  Output:
    - marked_duplicates.bam      : BAM with duplicates marked
    - duplicate_metrics.txt      : Picard MarkDuplicates metrics
    - duplicate_summary.txt      : parsed summary for final report
================================================================================
*/

process DUPLICATE_METRICS {
    label 'process_high'
    publishDir "${params.outdir}/duplicate_metrics/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_marked.bam"),
                     path("${meta.id}_marked.bai"),  emit: bam
    path "duplicate_metrics.txt",                     emit: metrics
    path "duplicate_summary.txt",                     emit: summary

    script:
    def prefix = "${meta.id}"
    """
    # ── Mark duplicates ───────────────────────────────────────────────────────
    picard MarkDuplicates \\
        I=${bam} \\
        O=${prefix}_marked.bam \\
        M=duplicate_metrics.txt \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT

    mv ${prefix}_marked.bai ${prefix}_marked.bai || true

    # ── Parse duplication rate ────────────────────────────────────────────────
    python3 - << 'PYEOF'
dup_rate = "NA"
est_library = "NA"
try:
    with open("duplicate_metrics.txt") as fh:
        lines = fh.readlines()
    for i, line in enumerate(lines):
        if line.startswith("ESTIMATED_LIBRARY_SIZE") or line.startswith("LIBRARY"):
            # Find the metrics header line
            pass
        if "PERCENT_DUPLICATION" in line:
            header = line.strip().split("\t")
            vals   = lines[i+1].strip().split("\t")
            idx_dup = header.index("PERCENT_DUPLICATION")
            dup_rate = float(vals[idx_dup])
            idx_lib = header.index("ESTIMATED_LIBRARY_SIZE") if "ESTIMATED_LIBRARY_SIZE" in header else -1
            if idx_lib >= 0:
                est_library = vals[idx_lib]
            break
except Exception as e:
    print(f"Warning: could not parse duplicate_metrics.txt: {e}")

threshold = ${params.max_duplication_rate}
flag = "FAIL" if isinstance(dup_rate, float) and dup_rate > threshold else "PASS"

with open("duplicate_summary.txt", "w") as out:
    out.write(f"step=duplicate_metrics\n")
    out.write(f"dataset=${meta.id}\n")
    out.write(f"max_duplication_rate={threshold}\n")
    out.write(f"duplication_rate={dup_rate}\n")
    out.write(f"estimated_library_size={est_library}\n")
    out.write(f"status={flag}\n")

print(f"Duplicate metrics: {dup_rate:.4f} ({float(dup_rate)*100:.1f}%) — {flag}")
PYEOF
    """
}
