/*
================================================================================
  MODULE: HETEROZYGOSITY
================================================================================
  Purpose:
    Calculate the observed heterozygosity rate for each sample and flag
    outliers. Outliers suggest DNA contamination (elevated heterozygosity)
    or inbreeding / sample degradation (reduced heterozygosity).

  Why this step:
    Heterozygosity is calculated as (N_observed_het / N_non_missing_genotypes).
    Samples deviating more than ±params.heterozygosity_sd standard deviations
    from the mean are excluded. Anderson et al. 2010 recommend ±3 SD as a
    standard threshold; tighten to ±2 SD for high-quality cohorts.

  Default threshold: params.heterozygosity_sd = 3

  How to change:
    nextflow run main.nf --heterozygosity_sd 2

  How to disable:
    nextflow run main.nf --run_heterozygosity false

  Important:
    Run on LD-pruned, common variants only to avoid confounding by regions of
    high LD or runs of homozygosity. This module performs LD pruning internally.

  Output:
    - heterozygosity.het         : PLINK heterozygosity table (all samples)
    - heterozygosity_outliers.txt : FID IID of outlier samples
    - heterozygosity_summary.txt  : mean, SD, cutoffs, removal count
    - heterozygosity_plot.png     : distribution plot (if R available)
================================================================================
*/

process HETEROZYGOSITY {
    label 'process_medium'
    publishDir "${params.outdir}/qc_tables", mode: params.publish_dir_mode, pattern: "*.{het,txt}"
    publishDir "${params.outdir}/qc_plots",  mode: params.publish_dir_mode, pattern: "*.png"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    path "heterozygosity_outliers.txt",  emit: outlier_samples
    path "heterozygosity_summary.txt",   emit: summary
    path "heterozygosity.het",           emit: het

    script:
    def prefix = "${meta.id}"
    """
    # ── LD pruning before heterozygosity calculation ──────────────────────────
    # LD-pruned, common SNPs give a cleaner estimate of the genome-wide het rate
    plink \\
        --bfile ${bed.baseName} \\
        --maf 0.05 \\
        --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \\
        --out prune_het \\
        --allow-no-sex

    # ── Calculate heterozygosity on pruned set ────────────────────────────────
    plink \\
        --bfile ${bed.baseName} \\
        --extract prune_het.prune.in \\
        --het \\
        --out heterozygosity \\
        --allow-no-sex

    # ── Identify outliers in shell (awk/python) ───────────────────────────────
    # .het columns: FID IID O(HOM) E(HOM) N(NM) F
    # Heterozygosity = (N(NM) - O(HOM)) / N(NM)
    python3 - << 'PYEOF'
import sys, math

data = []
with open("heterozygosity.het") as fh:
    next(fh)  # skip header
    for line in fh:
        parts = line.split()
        fid, iid = parts[0], parts[1]
        o_hom = float(parts[2])
        n_nm  = float(parts[4])
        if n_nm > 0:
            het_rate = (n_nm - o_hom) / n_nm
            data.append((fid, iid, het_rate))

if not data:
    sys.exit("ERROR: no data in heterozygosity.het")

rates = [d[2] for d in data]
mean_het = sum(rates) / len(rates)
sd_het   = math.sqrt(sum((r - mean_het)**2 for r in rates) / (len(rates) - 1))
sd_mult  = ${params.heterozygosity_sd}
lo = mean_het - sd_mult * sd_het
hi = mean_het + sd_mult * sd_het

outliers = [(fid, iid, r) for fid, iid, r in data if r < lo or r > hi]

with open("heterozygosity_outliers.txt", "w") as out:
    for fid, iid, r in outliers:
        out.write(f"{fid}\t{iid}\n")

with open("heterozygosity_summary.txt", "w") as out:
    out.write(f"step=heterozygosity\n")
    out.write(f"dataset=${meta.id}\n")
    out.write(f"sd_threshold=${params.heterozygosity_sd}\n")
    out.write(f"mean_het_rate={mean_het:.6f}\n")
    out.write(f"sd_het_rate={sd_het:.6f}\n")
    out.write(f"lower_cutoff={lo:.6f}\n")
    out.write(f"upper_cutoff={hi:.6f}\n")
    out.write(f"n_samples={len(data)}\n")
    out.write(f"n_outliers_removed={len(outliers)}\n")

print(f"Heterozygosity: {len(outliers)} outliers flagged "
      f"(mean={mean_het:.4f}, SD={sd_het:.4f}, "
      f"cutoffs=[{lo:.4f}, {hi:.4f}])")
PYEOF

    # ── Optional R plot ───────────────────────────────────────────────────────
    if command -v Rscript &>/dev/null; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
df <- read.table("heterozygosity.het", header=TRUE)
df\$het_rate <- (df\$N.NM. - df\$O.HOM.) / df\$N.NM.
mean_het <- mean(df\$het_rate, na.rm=TRUE)
sd_het   <- sd(df\$het_rate, na.rm=TRUE)
lo <- mean_het - ${params.heterozygosity_sd} * sd_het
hi <- mean_het + ${params.heterozygosity_sd} * sd_het
p <- ggplot(df, aes(x=het_rate)) +
    geom_histogram(bins=80, fill="steelblue", alpha=0.8) +
    geom_vline(xintercept=c(lo, hi), linetype="dashed", colour="red") +
    geom_vline(xintercept=mean_het, linetype="solid", colour="darkblue") +
    labs(title=paste0("Heterozygosity (±${params.heterozygosity_sd} SD cutoffs shown)"),
         x="Observed heterozygosity rate", y="Count") +
    theme_classic()
ggsave("heterozygosity_plot.png", p, width=8, height=5)
RSCRIPT
    fi
    """
}
