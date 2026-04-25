/*
================================================================================
  MODULE: FASTQC
================================================================================
  Purpose:
    Run FastQC on raw FASTQ files to check read quality before alignment.
    Generates per-base quality scores, adapter content, GC content, and
    duplication estimates. MultiQC aggregates results across samples.

  Why this step:
    Poor read quality, high adapter content, or unusual GC distributions
    indicate library preparation problems that downstream steps cannot fix.
    Detecting these early saves compute time and prevents misinterpreting
    downstream QC failures.

  How to disable:
    params.run_fastqc = false

  Output:
    - *.fastqc.html / *.fastqc.zip   : per-sample FastQC reports
    - multiqc_report.html            : aggregate report across samples
    - fastqc_summary.txt             : per-sample pass/warn/fail status
================================================================================
*/

process FASTQC {
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads_r1), path(reads_r2)   // reads_r2 may be []

    output:
    path "*.html",             emit: html
    path "*.zip",              emit: zip
    path "fastqc_summary.txt", emit: summary

    script:
    def r2_arg = reads_r2 ? "${reads_r2}" : ""
    def threads = task.cpus
    """
    # ── Run FastQC ────────────────────────────────────────────────────────────
    fastqc \\
        --threads ${threads} \\
        --outdir . \\
        ${reads_r1} ${r2_arg}

    # ── Parse FastQC summary for the report ───────────────────────────────────
    python3 - << 'PYEOF'
import glob, zipfile

results = {}
for zf in glob.glob("*.zip"):
    try:
        with zipfile.ZipFile(zf) as z:
            summary_name = [n for n in z.namelist() if n.endswith("summary.txt")]
            if summary_name:
                lines = z.read(summary_name[0]).decode()
                for line in lines.strip().split("\\n"):
                    status, test, fname = line.split("\t")
                    key = f"{fname}::{test}"
                    results[key] = status
    except Exception as e:
        print(f"Warning: could not parse {zf}: {e}")

n_fail = sum(1 for s in results.values() if s == "FAIL")
n_warn = sum(1 for s in results.values() if s == "WARN")
n_pass = sum(1 for s in results.values() if s == "PASS")

with open("fastqc_summary.txt", "w") as out:
    out.write("step=fastqc\\n")
    out.write(f"dataset=${meta.id}\\n")
    out.write(f"n_tests_pass={n_pass}\\n")
    out.write(f"n_tests_warn={n_warn}\\n")
    out.write(f"n_tests_fail={n_fail}\\n")
    for key, status in results.items():
        if status in ("FAIL", "WARN"):
            out.write(f"flag={status}::{key}\\n")

print(f"FastQC: {n_pass} PASS, {n_warn} WARN, {n_fail} FAIL across all tests")
PYEOF
    """
}
