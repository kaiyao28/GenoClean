/*
================================================================================
  MODULE: COVERAGE_QC
================================================================================
  Purpose:
    Calculate sequencing coverage metrics for WGS (genome-wide) or WES
    (on-target). Reports mean depth, fraction of bases at ≥10x, ≥20x, ≥30x.

  Why this step:
    Insufficient coverage is the most common reason for low variant call quality.
    GATK Best Practices recommend ≥20x mean depth for WGS and ≥30x for WES
    to achieve reliable variant calling. gnomAD uses coverage at multiple depth
    thresholds as a primary sample QC metric.

  Thresholds:
    params.min_mean_depth_wgs    = 20
    params.min_mean_depth_wes    = 30
    params.min_target_20x_fraction = 0.80  (fraction of targets covered ≥20x)

  How to change:
    nextflow run main.nf --min_mean_depth_wes 40

  How to disable:
    nextflow run main.nf --run_coverage_qc false

  Output:
    - coverage_summary_metrics.txt  : samtools depth / mosdepth metrics
    - coverage_summary.txt          : parsed summary for final report
    - coverage_plot.png             : depth distribution histogram (if R available)
================================================================================
*/

process COVERAGE_QC {
    label 'process_medium'
    publishDir "${params.outdir}/coverage_qc/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path reference_fasta
    path target_intervals   // may be [] for WGS

    output:
    path "coverage_summary.txt",       emit: summary
    path "coverage_metrics.txt",       emit: metrics

    script:
    def mode          = meta.mode ?: params.mode
    def min_depth     = mode == "wes" ? params.min_mean_depth_wes : params.min_mean_depth_wgs
    def has_intervals = target_intervals instanceof List ? false : target_intervals.toString() != "[]"
    def interval_arg  = has_intervals ? "--by ${target_intervals}" : ""
    """
    # ── Calculate coverage with mosdepth ──────────────────────────────────────
    # mosdepth is fast and outputs quantile coverage files
    mosdepth \\
        --threads ${task.cpus} \\
        --quantize 0:1:10:20:30:100: \\
        ${interval_arg} \\
        coverage \\
        ${bam}

    # ── Parse mosdepth output ─────────────────────────────────────────────────
    python3 - << 'PYEOF'
import gzip, os

mode         = "${mode}"
min_depth    = ${min_depth}
frac_20x_thr = ${params.min_target_20x_fraction}

# mosdepth summary file
mean_depth  = "NA"
pct_10x     = "NA"
pct_20x     = "NA"
pct_30x     = "NA"

try:
    with open("coverage.mosdepth.summary.txt") as fh:
        lines = fh.readlines()
    header = lines[0].strip().split("\t")
    # Find total or target region row
    for line in lines[1:]:
        parts = line.strip().split("\t")
        region = parts[0]
        if region in ("total", "total_region", "genome") or (
                "${has_intervals}" == "true" and "_region" in region):
            mean_col = header.index("mean") if "mean" in header else 3
            mean_depth = float(parts[mean_col])
            break
except Exception as e:
    print(f"Warning parsing mosdepth summary: {e}")

# Parse quantized regions to compute coverage fractions
# mosdepth quantize output: chr start end label
counts = {"0:1":0, "1:10":0, "10:20":0, "20:30":0, "30:100":0}
total_bases = 0
try:
    qfile = "coverage.quantized.bed.gz"
    if not os.path.exists(qfile):
        qfile = "coverage.quantized.bed"
    opener = gzip.open if qfile.endswith(".gz") else open
    with opener(qfile, "rt") as fh:
        for line in fh:
            chrom, start, end, label = line.strip().split("\t")
            bases = int(end) - int(start)
            if label in counts:
                counts[label] += bases
            total_bases += bases
except Exception as e:
    print(f"Warning parsing quantized bed: {e}")

if total_bases > 0:
    bases_0x  = counts.get("0:1", 0)
    bases_lt10 = counts.get("0:1",0) + counts.get("1:10",0)
    pct_10x   = 1 - (bases_lt10 / total_bases)
    bases_lt20 = bases_lt10 + counts.get("10:20",0)
    pct_20x   = 1 - (bases_lt20 / total_bases)
    bases_lt30 = bases_lt20 + counts.get("20:30",0)
    pct_30x   = 1 - (bases_lt30 / total_bases)

# Determine pass/fail
depth_ok  = isinstance(mean_depth, float) and mean_depth >= min_depth
frac_ok   = isinstance(pct_20x, float) and pct_20x >= frac_20x_thr
status    = "PASS" if (depth_ok and frac_ok) else "FAIL"

with open("coverage_summary.txt", "w") as out:
    out.write(f"step=coverage_qc\n")
    out.write(f"dataset=${meta.id}\n")
    out.write(f"mode={mode}\n")
    out.write(f"min_mean_depth_threshold={min_depth}\n")
    out.write(f"min_20x_fraction_threshold={frac_20x_thr}\n")
    out.write(f"mean_depth={mean_depth}\n")
    out.write(f"fraction_10x={pct_10x}\n")
    out.write(f"fraction_20x={pct_20x}\n")
    out.write(f"fraction_30x={pct_30x}\n")
    out.write(f"status={status}\n")

# Machine-readable metrics table
with open("coverage_metrics.txt", "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"mean_depth\t{mean_depth}\n")
    out.write(f"fraction_10x\t{pct_10x}\n")
    out.write(f"fraction_20x\t{pct_20x}\n")
    out.write(f"fraction_30x\t{pct_30x}\n")
    out.write(f"status\t{status}\n")

print(f"Coverage QC: mean depth={mean_depth}, 20x fraction={pct_20x} — {status}")
PYEOF

    # ── Optional depth histogram plot ─────────────────────────────────────────
    if command -v Rscript &>/dev/null && [ -f coverage.mosdepth.global.dist.txt ]; then
        Rscript - << 'RSCRIPT'
df <- read.table("coverage.mosdepth.global.dist.txt", header=FALSE,
                 col.names=c("region","depth","fraction"))
df <- df[df\$region == "total",]
library(ggplot2)
p <- ggplot(df[df\$depth <= 100,], aes(x=depth, y=fraction)) +
    geom_line(colour="steelblue") +
    geom_vline(xintercept=c(10,20,30), linetype="dashed", colour="red") +
    labs(title=paste0("Coverage distribution — ${meta.id}"),
         x="Depth", y="Fraction of bases at ≥ depth") +
    theme_classic()
ggsave("coverage_plot.png", p, width=8, height=5)
RSCRIPT
    fi
    """
}
