#!/usr/bin/env nextflow
/*
================================================================================
  SNP Array QC — inspect.nf
================================================================================
  Pre-QC inspection: computes all QC metric distributions WITHOUT applying
  any filters. Designed to run before main.nf so you can choose thresholds
  that are appropriate for your specific dataset.

  Outputs (in --outdir):
    params_template.yaml   — all default params pre-filled and annotated with
                             observed data stats; edit and pass to main.nf
    inspect_report.html    — browsable summary: plots + observed stats table
    qc_plots/              — individual PNG plots for each metric

  Usage:
    # Step 1 — inspect (no filtering)
    nextflow run snp_array_qc/inspect.nf \
      --bfile data/raw/genotypes \
      --outdir results/inspect

    # Step 2 — edit the template (optional; defaults work for most datasets)
    vi results/inspect/params_template.yaml

    # Step 3 — run full QC
    nextflow run snp_array_qc/main.nf \
      -params-file results/inspect/params_template.yaml

    # Or skip inspect entirely and run with all defaults in one line:
    nextflow run snp_array_qc/main.nf --bfile data/raw/genotypes
================================================================================
*/

nextflow.enable.dsl = 2

include { INPUT_CHECK } from './modules/input_check'

def validateParams() {
    if (!params.bfile) error "ERROR: --bfile is required."
    if (!file("${params.bfile}.bed").exists()) error "ERROR: ${params.bfile}.bed not found"
    if (!file("${params.bfile}.bim").exists()) error "ERROR: ${params.bfile}.bim not found"
    if (!file("${params.bfile}.fam").exists()) error "ERROR: ${params.bfile}.fam not found"
}

// ══════════════════════════════════════════════════════════════════════════════
//  WORKFLOW
// ══════════════════════════════════════════════════════════════════════════════
workflow {

    validateParams()

    log.info """
    ================================================================
    SNP Array QC — Inspection (no filtering)
    ================================================================
    PLINK prefix  : ${params.bfile}
    Output dir    : ${params.outdir}
    ================================================================
    No filters will be applied. Outputs:
      params_template.yaml  ← edit thresholds, then pass to main.nf
      inspect_report.html   ← distribution plots + observed stats
      qc_plots/*.png        ← individual metric plots
    ================================================================
    """.stripIndent()

    def prefix = file(params.bfile).name
    ch_plink = Channel.of([
        [ id: prefix ],
        file("${params.bfile}.bed"),
        file("${params.bfile}.bim"),
        file("${params.bfile}.fam")
    ])

    INPUT_CHECK(ch_plink)

    // Three stat processes run in parallel on the same input
    INSPECT_STATS(INPUT_CHECK.out.plink)
    INSPECT_RELATEDNESS(INPUT_CHECK.out.plink)
    INSPECT_PCA(INPUT_CHECK.out.plink)

    ch_all_plots = INSPECT_STATS.out.plots
        .mix(INSPECT_RELATEDNESS.out.plots)
        .mix(INSPECT_PCA.out.plots)

    PARAMS_TEMPLATE(
        INPUT_CHECK.out.summary,
        INSPECT_STATS.out.stats_summary,
        INSPECT_RELATEDNESS.out.ibd_summary,
        INSPECT_PCA.out.pca_summary,
        ch_all_plots.collect().ifEmpty([])
    )
}

workflow.onComplete {
    def status = workflow.success ? "COMPLETE" : "FAILED"
    def outdir = params.outdir
    log.info """
    ================================================================
    Inspection ${status} — no data was filtered
    ================================================================
    1. Open in your browser:
       ${outdir}/inspect_report.html

    2. Edit thresholds if needed (optional):
       ${outdir}/params_template.yaml
       — Each threshold is annotated with observed stats from your data.
       — Lines marked RECOMMENDED: fill in if you have the file.
       — Lines marked null: leave as-is if not needed.

    3. Run QC:
       nextflow run snp_array_qc/main.nf \\
         -params-file ${outdir}/params_template.yaml \\
         -profile <your_profile> -resume
    ================================================================
    """.stripIndent()
}

// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS: INSPECT_STATS
//  Per-sample: missingness, heterozygosity, sex check F-statistic.
//  Per-variant: missingness, MAF, HWE p-values.
//  All computed on unfiltered data; no --make-bed calls.
// ══════════════════════════════════════════════════════════════════════════════
process INSPECT_STATS {
    label 'process_high'
    publishDir "${params.outdir}/qc_plots", mode: params.publish_dir_mode, pattern: "*.png"
    publishDir "${params.outdir}/logs",     mode: params.publish_dir_mode, pattern: "inspect_stats.tsv"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    path "inspect_stats.tsv",            emit: stats_summary
    path "*.png",   optional: true,      emit: plots

    script:
    def prefix = "${meta.id}"
    """
    # ── Missingness (sample + variant) ───────────────────────────────────────
    plink --bfile ${bed.baseName} --missing --out ${prefix}_miss --allow-no-sex

    # ── MAF distribution ─────────────────────────────────────────────────────
    plink --bfile ${bed.baseName} --freq --out ${prefix}_freq --allow-no-sex

    # ── HWE p-values ─────────────────────────────────────────────────────────
    plink --bfile ${bed.baseName} --hardy --out ${prefix}_hwe --allow-no-sex

    # ── Sex check F-statistic ─────────────────────────────────────────────────
    plink --bfile ${bed.baseName} --check-sex --out ${prefix}_sex --allow-no-sex || true

    # ── Heterozygosity (LD-pruned, common SNPs only) ──────────────────────────
    plink \\
        --bfile ${bed.baseName} \\
        --maf 0.05 --geno 0.01 --hwe 0.001 \\
        --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \\
        --out prune_het --allow-no-sex
    plink \\
        --bfile ${bed.baseName} \\
        --extract prune_het.prune.in \\
        --het --out ${prefix}_het --allow-no-sex

    # ── Compute key percentiles (Python) ──────────────────────────────────────
    python3 - << 'PYEOF'
import math, os

def percentile(vals, p):
    if not vals: return 0.0
    idx = max(0, min(int(len(vals) * p / 100), len(vals) - 1))
    return vals[idx]

def read_col(fname, col, header=True):
    vals = []
    if not os.path.exists(fname): return vals
    with open(fname) as fh:
        if header: next(fh)
        for line in fh:
            parts = line.split()
            if len(parts) > col:
                try: vals.append(float(parts[col]))
                except ValueError: pass
    return vals

stats = {}

# Sample missingness
s_miss = sorted(read_col("${prefix}_miss.imiss", 4))
stats["sample_miss_median"] = f"{percentile(s_miss, 50):.6f}"
stats["sample_miss_p95"]    = f"{percentile(s_miss, 95):.6f}"
stats["sample_miss_max"]    = f"{s_miss[-1]:.6f}" if s_miss else "0.0"
stats["n_samples"]          = str(len(s_miss))

# Variant missingness
v_miss = sorted(read_col("${prefix}_miss.lmiss", 4))
stats["variant_miss_median"] = f"{percentile(v_miss, 50):.6f}"
stats["variant_miss_p95"]    = f"{percentile(v_miss, 95):.6f}"
stats["variant_miss_max"]    = f"{v_miss[-1]:.6f}" if v_miss else "0.0"
stats["n_variants"]          = str(len(v_miss))

# MAF
mafs = sorted(read_col("${prefix}_freq.frq", 4))
stats["maf_p5"]    = f"{percentile(mafs, 5):.6f}"
stats["maf_p25"]   = f"{percentile(mafs, 25):.6f}"
stats["maf_median"] = f"{percentile(mafs, 50):.6f}"

# HWE — separate autosomes and chrX
hwe_auto_p, hwe_chrx_p = [], []
has_chrx = False
hwe_file = "${prefix}_hwe.hwe"
if os.path.exists(hwe_file):
    with open(hwe_file) as fh:
        next(fh)
        for line in fh:
            parts = line.split()
            if len(parts) < 9: continue
            try:
                chrom = int(parts[0])
                p     = float(parts[8])
                if chrom == 23:
                    has_chrx = True
                    hwe_chrx_p.append(p)
                elif chrom < 23:
                    hwe_auto_p.append(p)
            except (ValueError, IndexError):
                pass
stats["hwe_auto_min_p"]       = f"{min(hwe_auto_p):.2e}" if hwe_auto_p else "1.0"
stats["hwe_auto_n_below_1e6"] = str(sum(1 for p in hwe_auto_p if p < 1e-6))
stats["hwe_chrx_min_p"]       = f"{min(hwe_chrx_p):.2e}" if hwe_chrx_p else "1.0"
stats["has_chrx"]             = "true" if has_chrx else "false"

# Heterozygosity
het_rates = []
het_file = "${prefix}_het.het"
if os.path.exists(het_file):
    with open(het_file) as fh:
        next(fh)
        for line in fh:
            p = line.split()
            if len(p) >= 5:
                try:
                    o_hom = float(p[2])
                    n     = float(p[4])
                    if n > 0: het_rates.append((n - o_hom) / n)
                except (ValueError, IndexError):
                    pass
if het_rates:
    mean = sum(het_rates) / len(het_rates)
    sd   = math.sqrt(sum((x - mean)**2 for x in het_rates) / max(len(het_rates) - 1, 1))
else:
    mean = sd = 0.0
stats["het_mean"]      = f"{mean:.6f}"
stats["het_sd"]        = f"{sd:.6f}"
stats["het_lower_3sd"] = f"{mean - 3*sd:.6f}"
stats["het_upper_3sd"] = f"{mean + 3*sd:.6f}"

# Sex check F-statistic
n_male = n_female = n_ambig = 0
sex_file = "${prefix}_sex.sexcheck"
if os.path.exists(sex_file):
    with open(sex_file) as fh:
        next(fh)
        for line in fh:
            parts = line.split()
            if len(parts) < 6: continue
            try:
                f = float(parts[5])
                if f < 0.2:        n_female += 1
                elif f > 0.8:      n_male += 1
                else:              n_ambig += 1
            except (ValueError, IndexError):
                pass
stats["sex_n_male"]      = str(n_male)
stats["sex_n_female"]    = str(n_female)
stats["sex_n_ambiguous"] = str(n_ambig)

# Case-control status from FAM
n_cases = n_controls = 0
with open("${fam}") as fh:
    for line in fh:
        p = line.split()
        if len(p) >= 6:
            if p[5] == "2":   n_cases += 1
            elif p[5] == "1": n_controls += 1
stats["n_cc_cases"]    = str(n_cases)
stats["n_cc_controls"] = str(n_controls)

with open("inspect_stats.tsv", "w") as out:
    for k, v in stats.items():
        out.write(f"{k}\\t{v}\\n")

print(f"Inspect stats: {len(stats)} metrics written")
PYEOF

    # ── Plots (R) ─────────────────────────────────────────────────────────────
    if command -v Rscript &>/dev/null; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
theme_set(theme_classic())
thr_col <- "red"

safe_read <- function(f, ...) {
    if (!file.exists(f)) return(NULL)
    tryCatch(read.table(f, header=TRUE, ...), error=function(e) NULL)
}

# Sample missingness
df <- safe_read("${prefix}_miss.imiss")
if (!is.null(df)) {
    p <- ggplot(df, aes(x=F_MISS)) +
        geom_histogram(bins=80, fill="steelblue", alpha=0.8) +
        geom_vline(xintercept=0.02, linetype="dashed", colour=thr_col) +
        annotate("text", x=0.022, y=Inf, label="default 0.02",
                 vjust=2, hjust=0, size=3, colour=thr_col) +
        labs(title="Sample missingness — pre-QC",
             subtitle="Dashed: default threshold (0.02). Samples to the right would be removed.",
             x="Proportion missing (F_MISS)", y="Count")
    ggsave("inspect_sample_miss.png", p, width=8, height=5)
}

# Variant missingness
df <- safe_read("${prefix}_miss.lmiss")
if (!is.null(df)) {
    p <- ggplot(df, aes(x=F_MISS)) +
        geom_histogram(bins=80, fill="steelblue", alpha=0.8) +
        geom_vline(xintercept=0.02, linetype="dashed", colour=thr_col) +
        annotate("text", x=0.022, y=Inf, label="default 0.02",
                 vjust=2, hjust=0, size=3, colour=thr_col) +
        labs(title="Variant missingness — pre-QC",
             subtitle="Dashed: default threshold (0.02). Variants to the right would be removed.",
             x="Proportion missing (F_MISS)", y="Count")
    ggsave("inspect_variant_miss.png", p, width=8, height=5)
}

# MAF
df <- safe_read("${prefix}_freq.frq")
if (!is.null(df)) {
    p <- ggplot(df, aes(x=MAF)) +
        geom_histogram(bins=100, fill="steelblue", alpha=0.8) +
        geom_vline(xintercept=${params.maf}, linetype="dashed", colour=thr_col) +
        annotate("text", x=${params.maf}+0.003, y=Inf,
                 label=paste0("default ", ${params.maf}),
                 vjust=2, hjust=0, size=3, colour=thr_col) +
        labs(title="Minor allele frequency — pre-QC",
             subtitle="Dashed: default MAF threshold. Variants to the left would be removed.",
             x="Minor allele frequency", y="Count") +
        scale_x_continuous(limits=c(0, 0.5))
    ggsave("inspect_maf.png", p, width=8, height=5)
}

# HWE autosomes
df <- safe_read("${prefix}_hwe.hwe")
if (!is.null(df) && "P" %in% names(df)) {
    df_a <- df[df$CHR < 23 & !is.na(df$P) & df$P > 0, ]
    if (nrow(df_a) > 0) {
        df_a$logp <- -log10(df_a$P)
        thr <- -log10(${params.hwe_p})
        p <- ggplot(df_a, aes(x=logp)) +
            geom_histogram(bins=80, fill="steelblue", alpha=0.8) +
            geom_vline(xintercept=thr, linetype="dashed", colour=thr_col) +
            annotate("text", x=thr+0.2, y=Inf,
                     label=paste0("default (p=", ${params.hwe_p}, ")"),
                     vjust=2, hjust=0, size=3, colour=thr_col) +
            labs(title="HWE -log10(p) — autosomes, pre-QC",
                 subtitle="Dashed: default threshold. Variants to the right would be removed.",
                 x="-log10(HWE p-value)", y="Count")
        ggsave("inspect_hwe.png", p, width=8, height=5)
    }
}

# Heterozygosity
df <- safe_read("${prefix}_het.het")
if (!is.null(df)) {
    df$HET_RATE <- (df$N.NM. - df$O.HOM.) / df$N.NM.
    mn  <- mean(df$HET_RATE, na.rm=TRUE)
    sdd <- sd(df$HET_RATE, na.rm=TRUE)
    lo  <- mn - 3*sdd
    hi  <- mn + 3*sdd
    p <- ggplot(df, aes(x=HET_RATE)) +
        geom_histogram(bins=80, fill="steelblue", alpha=0.8) +
        geom_vline(xintercept=c(lo, hi), linetype="dashed", colour=thr_col) +
        annotate("text", x=lo-0.001, y=Inf, label="-3 SD",
                 vjust=2, hjust=1, size=3, colour=thr_col) +
        annotate("text", x=hi+0.001, y=Inf, label="+3 SD",
                 vjust=2, hjust=0, size=3, colour=thr_col) +
        labs(title="Heterozygosity rate — pre-QC",
             subtitle="Dashed: ±3 SD (default). Samples outside the lines would be removed.",
             x="Heterozygosity rate", y="Count")
    ggsave("inspect_heterozygosity.png", p, width=8, height=5)
}

# Sex check F-statistic
df <- safe_read("${prefix}_sex.sexcheck")
if (!is.null(df) && "F" %in% names(df)) {
    df2 <- df[!is.na(df$F), ]
    if (nrow(df2) > 0) {
        lo_f <- ${params.sex_check_f_lower_female}
        hi_f <- ${params.sex_check_f_upper_male}
        p <- ggplot(df2, aes(x=F)) +
            geom_histogram(bins=80, fill="steelblue", alpha=0.8) +
            geom_vline(xintercept=c(lo_f, hi_f), linetype="dashed", colour=thr_col) +
            annotate("text", x=lo_f-0.01, y=Inf, label="female cutoff",
                     vjust=2, hjust=1, size=3, colour=thr_col) +
            annotate("text", x=hi_f+0.01, y=Inf, label="male cutoff",
                     vjust=2, hjust=0, size=3, colour=thr_col) +
            labs(title="X-chromosome F-statistic — pre-QC",
                 subtitle="Values between dashed lines are ambiguous. Left = female, right = male.",
                 x="F statistic", y="Count") +
            scale_x_continuous(limits=c(-0.2, 1.2))
        ggsave("inspect_sex_check.png", p, width=8, height=5)
    }
}
RSCRIPT
    fi
    """
}

// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS: INSPECT_RELATEDNESS
//  Lightweight IBD scan. Flags pairs above pi_hat 0.125 (3rd-degree threshold).
//  Shows how many pairs would be removed at various thresholds.
// ══════════════════════════════════════════════════════════════════════════════
process INSPECT_RELATEDNESS {
    label 'process_high'
    publishDir "${params.outdir}/qc_plots", mode: params.publish_dir_mode, pattern: "*.png"
    publishDir "${params.outdir}/logs",     mode: params.publish_dir_mode, pattern: "inspect_ibd_summary.tsv"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    path "inspect_ibd_summary.tsv",  emit: ibd_summary
    path "*.png",   optional: true,  emit: plots

    script:
    """
    plink \\
        --bfile ${bed.baseName} \\
        --maf 0.05 --geno 0.02 \\
        --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \\
        --out prune_ibd --allow-no-sex

    plink \\
        --bfile ${bed.baseName} \\
        --extract prune_ibd.prune.in \\
        --genome --min 0.125 \\
        --out ibd_scan \\
        --allow-no-sex

    python3 - << 'PYEOF'
import os

thresholds = [
    ("above_3rd_degree",  0.125),
    ("above_default",     0.1875),
    ("above_2nd_degree",  0.25),
    ("above_1st_degree",  0.5),
    ("duplicate_or_MZ",   0.9),
]
counts = {k: 0 for k, _ in thresholds}

genome_file = "ibd_scan.genome"
if os.path.exists(genome_file):
    with open(genome_file) as fh:
        next(fh)
        for line in fh:
            parts = line.split()
            if len(parts) < 10: continue
            try:
                pi = float(parts[9])
                for label, thr in thresholds:
                    if pi >= thr:
                        counts[label] += 1
            except (ValueError, IndexError):
                pass

with open("inspect_ibd_summary.tsv", "w") as out:
    for label, _ in thresholds:
        out.write(f"{label}\\t{counts[label]}\\n")

for label, _ in thresholds:
    print(f"  IBD pairs {label}: {counts[label]}")
PYEOF

    if command -v Rscript &>/dev/null && [ -f ibd_scan.genome ]; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
df <- tryCatch(read.table("ibd_scan.genome", header=TRUE), error=function(e) NULL)
if (!is.null(df) && "PI_HAT" %in% names(df) && nrow(df) > 0) {
    thrs <- data.frame(
        x     = c(0.125, 0.1875, 0.25, 0.5),
        label = c("3rd degree", "default (0.1875)", "2nd degree", "1st degree")
    )
    p <- ggplot(df, aes(x=PI_HAT)) +
        geom_histogram(bins=60, fill="steelblue", alpha=0.8) +
        geom_vline(data=thrs, aes(xintercept=x), linetype="dashed", colour="red", alpha=0.7) +
        geom_text(data=thrs, aes(x=x, label=label), y=Inf,
                  vjust=2, hjust=-0.05, size=2.8, colour="red", angle=90) +
        labs(title="Pairwise IBD (PI_HAT) — pre-QC, pairs > 0.125 only",
             subtitle="Default removal threshold: 0.1875. Pairs to the right of that line would be flagged.",
             x="PI_HAT", y="Count") +
        scale_x_continuous(limits=c(0.125, 1.0)) +
        theme_classic()
    ggsave("inspect_relatedness.png", p, width=8, height=5)
}
RSCRIPT
    fi
    """
}

// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS: INSPECT_PCA
//  Quick 10-PC PCA to reveal population structure before filtering.
// ══════════════════════════════════════════════════════════════════════════════
process INSPECT_PCA {
    label 'process_high'
    publishDir "${params.outdir}/qc_plots", mode: params.publish_dir_mode, pattern: "*.png"
    publishDir "${params.outdir}/logs",     mode: params.publish_dir_mode, pattern: "inspect_pca_summary.tsv"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    path "inspect_pca_summary.tsv", emit: pca_summary
    path "*.png",  optional: true,  emit: plots

    script:
    """
    plink \\
        --bfile ${bed.baseName} \\
        --maf 0.05 --geno 0.01 --hwe 0.001 \\
        --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \\
        --out prune_pca --allow-no-sex

    plink \\
        --bfile ${bed.baseName} \\
        --extract prune_pca.prune.in \\
        --pca 10 \\
        --out pca \\
        --allow-no-sex

    python3 - << 'PYEOF'
import os
pve = []
if os.path.exists("pca.eigenval"):
    with open("pca.eigenval") as fh:
        vals = [float(l.strip()) for l in fh if l.strip()]
    total = sum(vals) or 1.0
    pve = [v / total * 100 for v in vals]
with open("inspect_pca_summary.tsv", "w") as out:
    for i, v in enumerate(pve[:10]):
        out.write(f"pc{i+1}_pve\\t{v:.2f}\\n")
print(f"PCA complete: {len(pve)} PCs, PC1={pve[0]:.1f}% variance" if pve else "PCA complete")
PYEOF

    if command -v Rscript &>/dev/null && [ -f pca.eigenvec ]; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
pcs <- read.table("pca.eigenvec", header=FALSE)
colnames(pcs)[1:2] <- c("FID","IID")
colnames(pcs)[3:ncol(pcs)] <- paste0("PC", seq_len(ncol(pcs)-2))
p <- ggplot(pcs, aes(x=PC1, y=PC2)) +
    geom_point(alpha=0.5, size=0.9, colour="steelblue") +
    labs(title="PCA — PC1 vs PC2 (pre-QC, unfiltered)",
         subtitle="No reference panel. Add --reference_panel in main.nf for ancestry-labelled PCA.",
         x="PC1", y="PC2") +
    theme_classic()
ggsave("inspect_pca.png", p, width=8, height=6)
RSCRIPT
    fi
    """
}

// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS: PARAMS_TEMPLATE
//  Reads all inspection summaries. Writes:
//    params_template.yaml — annotated with observed stats from this dataset
//    inspect_report.html  — self-contained HTML with embedded plots + stats table
// ══════════════════════════════════════════════════════════════════════════════
process PARAMS_TEMPLATE {
    label 'process_low'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path input_summary
    path stats_summary
    path ibd_summary
    path pca_summary
    path plot_files

    output:
    path "params_template.yaml", emit: yaml
    path "inspect_report.html",  emit: report

    script:
    def bfile_abs   = file(params.bfile).toAbsolutePath().toString()
    def outdir_main = file(params.outdir).toAbsolutePath().parent.resolve("snp_array_qc").toString()
    // Pass all params as strings to preserve scientific notation formatting
    def p_smiss     = "${params.sample_missingness}"
    def p_vmiss     = "${params.variant_missingness}"
    def p_cc        = "${params.cc_miss_p}"
    def p_hwe       = "${params.hwe_p}"
    def p_hwex      = "${params.hwe_p_chrx}"
    def p_maf       = "${params.maf}"
    def p_het       = "${params.heterozygosity_sd}"
    def p_pi        = "${params.relatedness_pi_hat}"
    def p_sex_lo    = "${params.sex_check_f_lower_female}"
    def p_sex_hi    = "${params.sex_check_f_upper_male}"
    def p_pca_sd    = "${params.pca_outlier_sd}"
    def p_npcs      = "${params.n_pcs}"
    def p_ncov      = "${params.n_pcs_covariates}"
    def p_ldw       = "${params.ld_window}"
    def p_lds       = "${params.ld_step}"
    def p_ldr       = "${params.ld_r2}"
    def p_impr2     = "${params.imputation_r2}"
    """
    python3 - << 'PYEOF'
import os, base64, datetime

def read_kv(fname):
    d = {}
    if not os.path.exists(fname): return d
    with open(fname) as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            sep = "\\t" if "\\t" in line else "="
            parts = line.split(sep, 1)
            if len(parts) == 2:
                d[parts[0].strip()] = parts[1].strip()
    return d

def g(d, key, default="?"):
    return d.get(key, default)

inp  = read_kv("${input_summary}")
st   = read_kv("${stats_summary}")
ibd  = read_kv("${ibd_summary}")
pca  = read_kv("${pca_summary}")

n_samples  = g(inp, "n_samples")
n_variants = g(inp, "n_variants")
n_cases    = g(inp, "n_cases", "0")
n_controls = g(inp, "n_controls", "0")
chr_dist   = g(inp, "chromosome_distribution")
has_chrx   = g(st, "has_chrx", "false") == "true"

try:    is_cc = int(n_cases) > 0 and int(n_controls) > 0
except: is_cc = False

s_p95  = g(st, "sample_miss_p95",  "?")
s_max  = g(st, "sample_miss_max",  "?")
v_p95  = g(st, "variant_miss_p95", "?")
v_max  = g(st, "variant_miss_max", "?")
het_lo = g(st, "het_lower_3sd",    "?")
het_hi = g(st, "het_upper_3sd",    "?")
het_m  = g(st, "het_mean",         "?")
het_sd = g(st, "het_sd",           "?")
maf_p5  = g(st, "maf_p5",  "?")
maf_p25 = g(st, "maf_p25", "?")
maf_med = g(st, "maf_median", "?")
hwe_min = g(st, "hwe_auto_min_p", "?")
hwe_n6  = g(st, "hwe_auto_n_below_1e6", "?")
chrx_min = g(st, "hwe_chrx_min_p", "?")

n_3rd   = g(ibd, "above_3rd_degree",  "?")
n_def   = g(ibd, "above_default",     "?")
n_2nd   = g(ibd, "above_2nd_degree",  "?")
n_1st   = g(ibd, "above_1st_degree",  "?")
n_dup   = g(ibd, "duplicate_or_MZ",   "?")

pve_parts = [f"PC{i+1}={g(pca, f'pc{i+1}_pve', '?')}%" for i in range(5)]
pve_str   = ", ".join(pve_parts)

sex_m = g(st, "sex_n_male",      "?")
sex_f = g(st, "sex_n_female",    "?")
sex_a = g(st, "sex_n_ambiguous", "?")

today     = datetime.date.today().isoformat()
bfile_abs = "${bfile_abs}"
outdir_s  = "${outdir_main}"

# ── Threshold recommendations (annotate YAML) ─────────────────────────────────
try:
    sp95 = float(s_p95)
    if sp95 < 0.005:   smiss_note = f"  # observed 95th pct={s_p95}, max={s_max} — very clean; could tighten to 0.01"
    elif sp95 < 0.015: smiss_note = f"  # observed 95th pct={s_p95}, max={s_max} — default 0.02 looks appropriate"
    else:              smiss_note = f"  # observed 95th pct={s_p95}, max={s_max} — check inspect_sample_miss.png"
except:
    smiss_note = f"  # observed 95th pct={s_p95}, max={s_max}"

try:
    vp95 = float(v_p95)
    if vp95 < 0.005:   vmiss_note = f"  # observed 95th pct={v_p95}, max={v_max} — very clean; could tighten to 0.01"
    elif vp95 < 0.015: vmiss_note = f"  # observed 95th pct={v_p95}, max={v_max} — default 0.02 looks appropriate"
    else:              vmiss_note = f"  # observed 95th pct={v_p95}, max={v_max} — check inspect_variant_miss.png"
except:
    vmiss_note = f"  # observed 95th pct={v_p95}, max={v_max}"

cc_note   = f"  # {n_cases} cases, {n_controls} controls — HWE in controls only; cc_miss_p applied" if is_cc else "  # no case-control coding — cc_miss_p will not be applied"
chrx_note = f"  # chrX detected — separate permissive threshold applied" if has_chrx else "  # no chrX SNPs detected — this threshold will not be applied"

lines = [
    f"# SNP Array QC — params_template.yaml",
    f"# Generated: {today}",
    f"# Input:  {bfile_abs}",
    f"#         {n_samples} samples, {n_variants} variants",
    f"#         Chromosomes: {chr_dist}",
    f"#",
    f"# How to use:",
    f"#   1. Review the plots in: {os.path.dirname(bfile_abs)}/qc_plots/",
    f"#   2. Adjust thresholds below based on observed distributions",
    f"#   3. Run: nextflow run snp_array_qc/main.nf -params-file params_template.yaml",
    f"# ───────────────────────────────────────────────────────────────────────────",
    f"",
    f"# ── Input ──────────────────────────────────────────────────────────────────",
    f'bfile: "{bfile_abs}"',
    f'outdir: "{outdir_s}"',
    f"",
    f"# RECOMMENDED — fill in if you have these files:",
    f"reference_panel: null # 1000 Genomes prefix for ancestry-labelled PCA (e.g. /shared/data/1000G/1000G_hg19)",
    f"ld_regions: null      # high-LD regions file (MHC, chr8/17 inversions) in PLINK range format",
    f"",
    f"# OPTIONAL — leave as null if not needed:",
    f"pheno: null           # external phenotype file if your FAM column 6 is not in PLINK 1/2 coding",
    f"hapmap_info: null     # HapMap population info for colour-coded PCA plot",
    f"info_file: null       # imputation quality file — only for post-imputation QC",
    f"",
    f"# NOTE — FAM file column 6 phenotype coding:",
    f"#   1 = control, 2 = case, 0 or -9 = missing/no phenotype",
    f"#   If all values are 0/-9, HWE runs on all samples (not controls only) and",
    f"#   the differential missingness test (cc_miss_p) will not be applied.",
    f"#   Observed: {n_cases} cases, {n_controls} controls{' — case-control mode active' if is_cc else ' — cohort mode (no case/control coding detected)'}",
    f"",
    f"# ── Step switches ──────────────────────────────────────────────────────────",
    f"run_imputation_filter: false   # enable for post-imputation QC only",
    f"run_duplicate_check: true",
    f"run_sample_missingness: true",
    f"run_variant_missingness: true",
    f"run_sex_check: true",
    f"run_heterozygosity: true",
    f"run_relatedness: true",
    f"run_ancestry_pca: true",
    f"run_hwe: true",
    f"run_maf_filter: true",
    f"",
    f"# ── Thresholds ─────────────────────────────────────────────────────────────",
    f"sample_missingness: ${p_smiss}{smiss_note}",
    f"",
    f"variant_missingness: ${p_vmiss}{vmiss_note}",
    f"",
    f"cc_miss_p: ${p_cc}{cc_note}",
    f"",
    f"# HWE — observed autosomes: min p={hwe_min}, variants below 1e-6: {hwe_n6}",
    f"hwe_p: ${p_hwe}",
    f"hwe_p_chrx: ${p_hwex}{chrx_note}",
    f"",
    f"# MAF — observed: 5th pct={maf_p5}, 25th pct={maf_p25}, median={maf_med}",
    f"# For rare-variant analysis: set run_maf_filter: false or use 0.001",
    f"maf: ${p_maf}",
    f"",
    f"# Heterozygosity — observed mean={het_m}, SD={het_sd}",
    f"#   ±3 SD bounds: [{het_lo}, {het_hi}]",
    f"heterozygosity_sd: ${p_het}",
    f"",
    f"# Relatedness — IBD pairs above each threshold:",
    f"#   > 0.125 (3rd degree): {n_3rd}",
    f"#   > 0.1875 (default):   {n_def}",
    f"#   > 0.25 (2nd degree):  {n_2nd}",
    f"#   > 0.5 (1st degree):   {n_1st}",
    f"#   duplicates/MZ twins:  {n_dup}",
    f"relatedness_pi_hat: ${p_pi}",
    f"",
    f"# Sex check — observed: {sex_m} male, {sex_f} female, {sex_a} ambiguous (F between 0.2–0.8)",
    f"sex_check_f_lower_female: ${p_sex_lo}",
    f"sex_check_f_upper_male: ${p_sex_hi}",
    f"",
    f"# PCA — variance explained: {pve_str}",
    f"# Use the scree plot from main.nf to choose n_pcs_covariates",
    f"pca_outlier_sd: ${p_pca_sd}",
    f"n_pcs: ${p_npcs}",
    f"n_pcs_covariates: ${p_ncov}",
    f"",
    f"# ── LD pruning ─────────────────────────────────────────────────────────────",
    f"ld_window: ${p_ldw}",
    f"ld_step: ${p_lds}",
    f"ld_r2: ${p_ldr}",
    f"",
    f"# ── Imputation filter (only when run_imputation_filter: true) ──────────────",
    f"imputation_r2: ${p_impr2}",
    f"",
    f"# ── Output ─────────────────────────────────────────────────────────────────",
    f"publish_dir_mode: copy",
    f"keep_intermediate: false",
]

with open("params_template.yaml", "w") as out:
    out.write("\\n".join(lines) + "\\n")

# ── HTML report ───────────────────────────────────────────────────────────────
def embed_png(path):
    if not os.path.exists(path): return ""
    with open(path, "rb") as f:
        enc = base64.b64encode(f.read()).decode()
    return (f'<figure style="margin:1.5rem 0">'
            f'<img src="data:image/png;base64,{enc}" style="max-width:100%;border:1px solid #ddd">'
            f'<figcaption style="font-size:0.85rem;color:#666;margin-top:0.3rem">{os.path.basename(path)}</figcaption>'
            f'</figure>')

plots_html = ""
for fname in sorted(f for f in os.listdir(".") if f.endswith(".png")):
    plots_html += embed_png(fname)

rows = [
    ("Samples",                  n_samples),
    ("Variants",                 n_variants),
    ("Cases / Controls",         f"{n_cases} / {n_controls}" if is_cc else "not case-control"),
    ("Chromosomes",              chr_dist),
    ("chrX SNPs present",        "yes" if has_chrx else "no"),
    ("Sample miss — 95th pct",   s_p95),
    ("Sample miss — max",        s_max),
    ("Variant miss — 95th pct",  v_p95),
    ("Variant miss — max",       v_max),
    ("Heterozygosity mean ± SD", f"{het_m} ± {het_sd}"),
    ("Het ±3 SD bounds",         f"[{het_lo}, {het_hi}]"),
    ("HWE min p (autosomes)",    hwe_min),
    ("HWE variants below 1e-6",  hwe_n6),
    ("IBD pairs > 0.1875 (default)", n_def),
    ("IBD duplicates / MZ twins",    n_dup),
    ("Sex — male / female / ambig",  f"{sex_m} / {sex_f} / {sex_a}"),
    ("PCA variance PC1–PC5",     pve_str),
]
rows_html = "".join(f"<tr><td>{l}</td><td><b>{v}</b></td></tr>" for l, v in rows)

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SNP Array QC — Pre-QC Inspection</title>
<style>
  body {{ font-family: sans-serif; max-width: 1100px; margin: 2rem auto; padding: 0 1rem; color: #222; }}
  h1 {{ border-bottom: 2px solid #2c6fad; padding-bottom: 0.4rem; }}
  h2 {{ color: #2c6fad; margin-top: 2rem; }}
  table {{ border-collapse: collapse; width: 60%; margin-bottom: 1.5rem; }}
  th, td {{ border: 1px solid #ccc; padding: 0.4rem 0.7rem; text-align: left; }}
  th {{ background: #e8f0fa; }}
  .note {{ background: #fffbe6; border-left: 4px solid #f0a500; padding: 0.8rem 1rem; margin: 1rem 0; }}
  code {{ background: #f4f4f4; padding: 0.1rem 0.4rem; border-radius: 3px; font-size: 0.92em; }}
</style>
</head>
<body>
<h1>SNP Array QC — Pre-QC Inspection Report</h1>
<p>Generated: {today} &nbsp;|&nbsp; Input: <code>{bfile_abs}</code></p>

<div class="note">
  <b>No filters were applied.</b> Plots show raw data distributions with default thresholds marked.
  Edit <code>params_template.yaml</code> based on these distributions, then run:<br>
  <code>nextflow run snp_array_qc/main.nf -params-file params_template.yaml</code>
</div>

<h2>Dataset summary</h2>
<table>
<tr><th>Metric</th><th>Observed</th></tr>
{rows_html}
</table>

<h2>Distribution plots</h2>
<p>Red dashed lines show default thresholds. Adjust in <code>params_template.yaml</code>.</p>
{plots_html}
</body>
</html>"""

with open("inspect_report.html", "w") as out:
    out.write(html)

print("params_template.yaml written")
print("inspect_report.html written")
PYEOF
    """
}
