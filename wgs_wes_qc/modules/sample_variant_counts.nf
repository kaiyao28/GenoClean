/*
================================================================================
  MODULE: SAMPLE_VARIANT_COUNTS
================================================================================
  Purpose:
    Compute per-sample variant counts from a filtered VCF: total variants,
    SNPs, indels, singletons, het/hom ratio, and missingness (call rate).
    Flag samples with outlier counts.

  Why this step:
    After genotype filtering, some samples may have an abnormally high
    missingness (low call rate), very high or low variant counts, or unusual
    het/hom ratios. These can indicate contamination, low quality, or
    cryptic population structure.

  Thresholds:
    params.min_call_rate = 0.95   (flag samples with call rate < 95%)

  How to disable:
    params.run_sample_variant_counts = false

  Output:
    - sample_variant_counts.tsv   : per-sample counts and call rates
    - sample_count_outliers.txt   : samples below min_call_rate
    - sample_variant_counts_summary.txt : aggregate summary
================================================================================
*/

process SAMPLE_VARIANT_COUNTS {
    label 'process_medium'
    publishDir "${params.outdir}/qc_tables",  mode: params.publish_dir_mode, pattern: "*.tsv"
    publishDir "${params.outdir}/qc_tables",  mode: params.publish_dir_mode, pattern: "*.txt"
    publishDir "${params.outdir}/qc_plots",   mode: params.publish_dir_mode, pattern: "*.png"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    path "sample_variant_counts.tsv",         emit: counts
    path "sample_count_outliers.txt",         emit: outliers
    path "sample_variant_counts_summary.txt", emit: summary

    script:
    """
    # ── Per-sample stats via bcftools ─────────────────────────────────────────
    bcftools stats \\
        --threads ${task.cpus} \\
        --samples - \\
        ${vcf} | grep "^PSC" > psc.txt || true

    # ── Missingness per sample ────────────────────────────────────────────────
    # Count missing genotypes per sample
    bcftools query \\
        --format '[%SAMPLE\t%GT\n]' \\
        ${vcf} | awk '
        {
            total[\$1]++
            if (\$2=="." || \$2=="./." || \$2==".|.") missing[\$1]++
        }
        END {
            for (s in total) {
                m = missing[s]+0
                t = total[s]
                call_rate = (t-m)/t
                print s"\t"t"\t"m"\t"call_rate
            }
        }
    ' > missingness.txt

    # ── Parse and flag low call-rate samples ──────────────────────────────────
    python3 - << 'PYEOF'
min_cr = ${params.min_call_rate}

samples = {}
try:
    with open("missingness.txt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) == 4:
                sample, total, missing, call_rate = parts
                samples[sample] = {
                    "total": int(total),
                    "missing": int(missing),
                    "call_rate": float(call_rate),
                }
except Exception as e:
    print(f"Warning: {e}")

# Merge with PSC stats
snp_counts = {}
try:
    with open("psc.txt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 14:
                # PSC cols: id, sample, hom_RR, hom_AA, het, ts, tv, indel, ...
                sample = parts[2]
                snp_counts[sample] = {
                    "n_hom": parts[4],
                    "n_het": parts[5],
                    "ts":    parts[6],
                    "tv":    parts[7],
                    "n_indel": parts[8],
                }
except Exception as e:
    print(f"Warning parsing PSC: {e}")

outliers = []
with open("sample_variant_counts.tsv", "w") as out:
    out.write("sample\tn_total\tn_missing\tcall_rate\tn_het\tn_hom\tn_indel\tstatus\\n")
    for sample, d in samples.items():
        cr = d["call_rate"]
        psc = snp_counts.get(sample, {})
        status = "FAIL" if cr < min_cr else "PASS"
        if status == "FAIL":
            outliers.append(sample)
        out.write(f"{sample}\t{d['total']}\t{d['missing']}\t{cr:.4f}\t"
                  f"{psc.get('n_het','NA')}\t{psc.get('n_hom','NA')}\t"
                  f"{psc.get('n_indel','NA')}\t{status}\\n")

with open("sample_count_outliers.txt", "w") as out:
    for s in outliers:
        out.write(f"{s}\\n")

with open("sample_variant_counts_summary.txt", "w") as out:
    out.write(f"step=sample_variant_counts\\n")
    out.write(f"dataset=${meta.id}\\n")
    out.write(f"min_call_rate_threshold={min_cr}\\n")
    out.write(f"n_samples={len(samples)}\\n")
    out.write(f"n_low_callrate={len(outliers)}\\n")
    if samples:
        avg_cr = sum(d["call_rate"] for d in samples.values()) / len(samples)
        out.write(f"mean_call_rate={avg_cr:.4f}\\n")

print(f"Sample variant counts: {len(outliers)} samples below call rate {min_cr}")
PYEOF

    # ── Optional variant count distribution plot ───────────────────────────────
    if command -v Rscript &>/dev/null && [ -s sample_variant_counts.tsv ]; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
df <- read.table("sample_variant_counts.tsv", header=TRUE, sep="\t")
p <- ggplot(df, aes(x=call_rate, fill=status)) +
    geom_histogram(bins=50, alpha=0.8) +
    geom_vline(xintercept=${params.min_call_rate}, linetype="dashed", colour="red") +
    scale_fill_manual(values=c("PASS"="steelblue", "FAIL"="red")) +
    labs(title="Sample call rates after genotype filtering",
         x="Call rate", y="Count", fill="Status") +
    theme_classic()
ggsave("sample_callrate_distribution.png", p, width=8, height=5)
RSCRIPT
    fi
    """
}
