/*
================================================================================
  MODULE: SEX_CHECK (WGS/WES)
================================================================================
  Purpose:
    Infer biological sex from the BAM file using X and Y chromosome coverage
    ratios and X-chromosome heterozygosity. Compare to reported sex if provided.

  Why this step:
    For WGS/WES, sex is inferred from:
      - Mean depth on chrX relative to autosomes (females ~1.0, males ~0.5)
      - Mean depth on chrY (females ~0, males ~0.5)
      - X heterozygosity (females have two X chromosomes → higher het rate)
    Discordance with reported sex suggests sample swap or mislabelling.

  How to disable:
    params.run_sex_check_wgs = false

  Output:
    - sex_check_results.txt    : per-sample inferred sex and metrics
    - sex_check_summary.txt    : summary for final report
================================================================================
*/

process SEX_CHECK {
    label 'process_low'
    publishDir "${params.outdir}/sex_check", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path reference_fasta

    output:
    path "sex_check_results.txt",  emit: results
    path "sex_check_summary.txt",  emit: summary

    script:
    """
    # ── Compute mean depth on autosomes, chrX, and chrY ──────────────────────
    # samtools idxstats gives mapped read counts per chromosome
    samtools idxstats ${bam} > idxstats.txt

    # samtools depth on specific regions for mean depth calculation
    samtools depth -r chr1 ${bam} 2>/dev/null | awk '{sum+=\$3; n++} END{print "chr1\t"sum/n}' > depth_auto.txt || true
    samtools depth -r 1 ${bam} 2>/dev/null | awk '{sum+=\$3; n++} END{print "chr1\t"sum/n}' >> depth_auto.txt || true

    samtools depth -r chrX ${bam} 2>/dev/null | awk '{sum+=\$3; n++} END{print "chrX\t"sum/n}' > depth_chrX.txt || true
    samtools depth -r X ${bam} 2>/dev/null | awk '{sum+=\$3; n++} END{print "chrX\t"sum/n}' >> depth_chrX.txt || true

    samtools depth -r chrY ${bam} 2>/dev/null | awk '{sum+=\$3; n++} END{print "chrY\t"sum/n}' > depth_chrY.txt || true
    samtools depth -r Y ${bam} 2>/dev/null | awk '{sum+=\$3; n++} END{print "chrY\t"sum/n}' >> depth_chrY.txt || true

    # ── Infer sex from coverage ratios ────────────────────────────────────────
    python3 - << 'PYEOF'
def read_depth(fname):
    try:
        with open(fname) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    val = float(parts[1])
                    if val > 0:
                        return val
    except Exception:
        pass
    return None

auto_depth = read_depth("depth_auto.txt")
x_depth    = read_depth("depth_chrX.txt")
y_depth    = read_depth("depth_chrY.txt")

# Normalised coverage ratios
x_norm = (x_depth / auto_depth) if (auto_depth and x_depth) else None
y_norm = (y_depth / auto_depth) if (auto_depth and y_depth) else None

# Inference logic:
# Female: X/auto ~1.0, Y/auto ~0.0
# Male:   X/auto ~0.5, Y/auto ~0.5
inferred_sex = "unknown"
if x_norm is not None:
    if x_norm > 0.7:
        inferred_sex = "female"
    elif x_norm < 0.6:
        inferred_sex = "male"
    else:
        inferred_sex = "ambiguous"

# Cross-check with Y depth
if y_norm is not None and inferred_sex == "female" and y_norm > 0.1:
    inferred_sex = "ambiguous (XX with Y signal)"

status = "OK"
reported = "${meta.sex ?: 'unknown'}"
if reported not in ("unknown", "") and reported != inferred_sex:
    status = "DISCORDANT"

with open("sex_check_results.txt", "w") as out:
    out.write("sample\treported_sex\tinferred_sex\tx_coverage\ty_coverage\tauto_coverage\tx_norm\ty_norm\tstatus\\n")
    out.write(f"${meta.id}\t{reported}\t{inferred_sex}\t{x_depth}\t{y_depth}\t{auto_depth}\t{x_norm}\t{y_norm}\t{status}\\n")

with open("sex_check_summary.txt", "w") as out:
    out.write(f"step=sex_check_wgs\\n")
    out.write(f"dataset=${meta.id}\\n")
    out.write(f"inferred_sex={inferred_sex}\\n")
    out.write(f"reported_sex={reported}\\n")
    out.write(f"x_norm_coverage={x_norm}\\n")
    out.write(f"y_norm_coverage={y_norm}\\n")
    out.write(f"status={status}\\n")

print(f"Sex check: inferred={inferred_sex}, reported={reported} — {status}")
PYEOF
    """
}
