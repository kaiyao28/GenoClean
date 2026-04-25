/*
================================================================================
  MODULE: ALIGNMENT_METRICS
================================================================================
  Purpose:
    Compute alignment statistics from a BAM/CRAM file using samtools and Picard.
    Checks mapping rate, properly paired reads, insert size distribution, and
    chimeric read rate.

  Why this step:
    Low mapping rate, high chimeric read rate, or narrow insert size suggest
    alignment problems, incorrect reference genome, or library preparation issues.
    These metrics are independent of variant calling and should be checked before
    investing compute in downstream analysis.

  How to disable:
    params.run_alignment_metrics = false

  Output:
    - flagstat.txt                    : samtools flagstat
    - alignment_summary_metrics.txt   : Picard CollectAlignmentSummaryMetrics
    - insert_size_metrics.txt         : Picard CollectInsertSizeMetrics
    - insert_size_histogram.pdf       : Picard insert size histogram
    - alignment_summary.txt           : parsed summary for final report
================================================================================
*/

process ALIGNMENT_METRICS {
    label 'process_medium'
    publishDir "${params.outdir}/alignment_metrics/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path reference_fasta

    output:
    path "flagstat.txt",                 emit: flagstat
    path "alignment_summary_metrics.txt", emit: picard_alignment
    path "insert_size_metrics.txt",      emit: insert_size
    path "alignment_summary.txt",        emit: summary

    script:
    """
    # ── samtools flagstat ──────────────────────────────────────────────────────
    samtools flagstat --threads ${task.cpus} ${bam} > flagstat.txt

    # ── Picard alignment summary metrics ─────────────────────────────────────
    picard CollectAlignmentSummaryMetrics \\
        R=${reference_fasta} \\
        I=${bam} \\
        O=alignment_summary_metrics.txt

    # ── Picard insert size metrics ────────────────────────────────────────────
    picard CollectInsertSizeMetrics \\
        I=${bam} \\
        O=insert_size_metrics.txt \\
        H=insert_size_histogram.pdf \\
        W=600

    # ── Parse key metrics for report ──────────────────────────────────────────
    python3 - << 'PYEOF'
import re

def grep_flagstat(key, text):
    for line in text.split("\\n"):
        if key in line:
            m = re.search(r"([0-9]+)", line)
            return m.group(1) if m else "NA"
    return "NA"

with open("flagstat.txt") as fh:
    flag_text = fh.read()

total_reads   = grep_flagstat("in total", flag_text)
mapped_reads  = grep_flagstat("mapped (", flag_text)
properly_pair = grep_flagstat("properly paired", flag_text)
duplicates    = grep_flagstat("duplicates", flag_text)

# Extract mapping rate percentage
map_rate = "NA"
for line in flag_text.split("\\n"):
    if "mapped (" in line:
        m = re.search(r"[(]([0-9.]+)%[)]", line)
        if m:
            map_rate = m.group(1)
        break

# Extract median insert size from Picard output
median_insert = "NA"
try:
    with open("insert_size_metrics.txt") as fh:
        lines = fh.readlines()
    for i, line in enumerate(lines):
        if line.startswith("MEDIAN_INSERT_SIZE"):
            vals = lines[i+1].split("\t")
            header = line.split("\t")
            idx = header.index("MEDIAN_INSERT_SIZE")
            median_insert = vals[idx]
            break
except Exception:
    pass

with open("alignment_summary.txt", "w") as out:
    out.write(f"step=alignment_metrics\\n")
    out.write(f"dataset=${meta.id}\\n")
    out.write(f"total_reads={total_reads}\\n")
    out.write(f"mapped_reads={mapped_reads}\\n")
    out.write(f"mapping_rate_pct={map_rate}\\n")
    out.write(f"properly_paired={properly_pair}\\n")
    out.write(f"duplicates={duplicates}\\n")
    out.write(f"median_insert_size={median_insert}\\n")

print(f"Alignment metrics: {total_reads} total reads, {map_rate}% mapped, "
      f"median insert size {median_insert} bp")
PYEOF
    """
}
