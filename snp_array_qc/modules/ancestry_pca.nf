/*
================================================================================
  MODULE: ANCESTRY_PCA (SNP Array)
================================================================================
  Purpose:
    Perform principal component analysis (PCA) on the genotype data to detect
    population structure and identify ancestry outliers. Samples deviating
    >params.pca_outlier_sd SD from the population centre on any PC are flagged.

  Design distinction from WGS/WES PCA:
    This module performs LD pruning and PCA in a single step, WITHOUT a
    preceding PCA_VARIANT_SELECTION module. This is correct for SNP arrays
    because:
      - Arrays genotype common variants by design (usually MAF > 1–5%)
      - After standard QC (missingness, HWE, MAF), remaining variants are
        already suitable for PCA
      - No rare variants, indels, or coverage biases to exclude

    The WGS/WES pipeline (wgs_wes_qc/modules/ancestry_pca.nf) runs a
    PCA_VARIANT_SELECTION step BEFORE PCA to handle rare variants, indels,
    and uneven coverage from sequencing data. Do NOT add that step here.

  Why this step:
    Population stratification is a major confounder in GWAS. PCA identifies
    individuals whose ancestry deviates from the study population. If a reference
    panel (e.g. 1000 Genomes) is provided, PCs can be interpreted in the context
    of known ancestries and projected PCs produced.

  Default threshold: params.pca_outlier_sd = 6

  How to change:
    nextflow run main.nf --pca_outlier_sd 4

  How to disable:
    nextflow run main.nf --run_ancestry_pca false

  Output:
    - pca.eigenvec           : sample PC scores
    - pca.eigenval           : variance explained per PC
    - ancestry_outliers.txt  : FID IID of outlier samples
    - pca_summary.txt        : variance explained and outlier counts
    - pca_plot.png           : PC1 vs PC2 scatter plot (if R available)
================================================================================
*/

process ANCESTRY_PCA {
    label 'process_high'
    publishDir "${params.outdir}/qc_tables", mode: params.publish_dir_mode, pattern: "*.{eigenvec,eigenval,txt}"
    publishDir "${params.outdir}/qc_plots",  mode: params.publish_dir_mode, pattern: "*.png"
    publishDir "${params.outdir}",           mode: params.publish_dir_mode, pattern: "pca_covariates.txt"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val   reference_panel   // path prefix or empty list
    path  ld_regions        // optional: high-LD region file for --exclude range
    path  hapmap_info       // optional: hapmap_sample_info.txt for population labels

    output:
    path "ancestry_outliers.txt",  emit: outlier_samples
    path "pca_summary.txt",        emit: summary
    path "pca.eigenvec",           emit: eigenvec
    path "pca.eigenval",           emit: eigenval
    path "pca_covariates.txt",     emit: covariates
    path "*.png",                  optional: true, emit: plots

    script:
    def prefix      = "${meta.id}"
    def has_ref     = reference_panel instanceof List ? false : (reference_panel as String).length() > 0
    def ld_flag     = ld_regions instanceof List ? "" : "--exclude range ${ld_regions}"
    def hapmap_path = hapmap_info instanceof List ? "" : hapmap_info.toString()
    """
    # ── LD pruning ────────────────────────────────────────────────────────────
    # Exclude high-LD regions (MHC, inversions) before pruning so they don't
    # dominate the top PCs and obscure true ancestry variation
    plink \\
        --bfile ${bed.baseName} \\
        --maf 0.05 \\
        --geno 0.01 \\
        --hwe 0.001 \\
        ${ld_flag} \\
        --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \\
        --out prune_pca \\
        --allow-no-sex

    # ── PCA (study data only, or merged with reference) ───────────────────────
    ${has_ref ? """
    # Merge study data with reference panel
    # Align strand and allele coding before merging
    awk 'NR==FNR{a[\$2]=1; next} a[\$2]' prune_pca.prune.in ${bim} | \\
        awk '{print \$2}' > common_snps.txt

    plink \\
        --bfile ${reference_panel} \\
        --extract common_snps.txt \\
        --make-bed \\
        --out ref_subset \\
        --allow-no-sex

    plink \\
        --bfile ${bed.baseName} \\
        --extract prune_pca.prune.in \\
        --bmerge ref_subset \\
        --merge-mode 6 \\
        --make-bed \\
        --out merged_study_ref \\
        --allow-no-sex || true

    # If merge fails due to strand flips, flip and retry
    if [ -f merged_study_ref-merge.missnp ]; then
        plink --bfile ref_subset --flip merged_study_ref-merge.missnp \\
              --make-bed --out ref_flipped --allow-no-sex
        plink --bfile ${bed.baseName} --extract prune_pca.prune.in \\
              --bmerge ref_flipped --merge-mode 6 --make-bed \\
              --out merged_study_ref --allow-no-sex
    fi

    plink \\
        --bfile merged_study_ref \\
        --pca ${params.n_pcs} \\
        --out pca \\
        --allow-no-sex
    """ : """
    # No reference panel — PCA on study data only
    plink \\
        --bfile ${bed.baseName} \\
        --extract prune_pca.prune.in \\
        --pca ${params.n_pcs} \\
        --out pca \\
        --allow-no-sex
    """}

    # ── Identify ancestry outliers + write analysis-ready covariate file ────────
    python3 - << 'PYEOF'
import math, sys

# ── Read eigenvec (all samples, including reference panel if merged) ───────────
data = []
with open("pca.eigenvec") as fh:
    for line in fh:
        parts = line.split()
        fid, iid = parts[0], parts[1]
        pcs = list(map(float, parts[2:]))
        data.append((fid, iid, pcs))

n_pcs_computed = len(data[0][2]) if data else 0

# ── Outlier detection (run on all samples so reference panel anchors the mean) ─
means = [sum(d[2][i] for d in data) / len(data) for i in range(n_pcs_computed)]
sds   = [math.sqrt(sum((d[2][i] - means[i])**2 for d in data) / max(len(data) - 1, 1))
         for i in range(n_pcs_computed)]

sd_thr   = ${params.pca_outlier_sd}
outliers = []
for fid, iid, pcs in data:
    for i, (pc, m, s) in enumerate(zip(pcs, means, sds)):
        if s > 0 and abs(pc - m) > sd_thr * s:
            outliers.append((fid, iid))
            break

with open("ancestry_outliers.txt", "w") as out:
    for fid, iid in outliers:
        out.write(f"{fid}\t{iid}\\n")

# ── Variance explained ────────────────────────────────────────────────────────
eigenvals = []
try:
    with open("pca.eigenval") as fh:
        eigenvals = [float(l.strip()) for l in fh if l.strip()]
    total_var = sum(eigenvals)
    pve = [v / total_var * 100 for v in eigenvals]
except Exception:
    pve = []

# ── Covariate file: study samples only, first n_pcs_covariates PCs ───────────
# Study samples are those present in the input FAM file.
# When a reference panel is merged, the eigenvec contains all samples;
# filtering to the FAM ensures only study-cohort PCs go into the covariate file.
study_ids = set()
with open("${fam}") as fh:
    for line in fh:
        p = line.split()
        if len(p) >= 2:
            study_ids.add((p[0], p[1]))

n_cov = min(${params.n_pcs_covariates}, n_pcs_computed)
if ${params.n_pcs_covariates} > n_pcs_computed:
    print(f"WARNING: n_pcs_covariates (${params.n_pcs_covariates}) > PCs computed "
          f"({n_pcs_computed}); covariate file will contain {n_cov} PCs")

header = "FID\\tIID\\t" + "\\t".join(f"PC{i+1}" for i in range(n_cov))
n_study_written = 0
with open("pca_covariates.txt", "w") as out:
    out.write(header + "\\n")
    for fid, iid, pcs in data:
        if (fid, iid) in study_ids:
            vals = "\\t".join(f"{pcs[i]:.8f}" for i in range(n_cov))
            out.write(f"{fid}\\t{iid}\\t{vals}\\n")
            n_study_written += 1

with open("pca_summary.txt", "w") as out:
    out.write(f"step=ancestry_pca\\n")
    out.write(f"dataset=${meta.id}\\n")
    out.write(f"pca_outlier_sd=${params.pca_outlier_sd}\\n")
    out.write(f"n_pcs_computed={n_pcs_computed}\\n")
    out.write(f"n_pcs_covariates={n_cov}\\n")
    out.write(f"n_samples_total={len(data)}\\n")
    out.write(f"n_samples_study={n_study_written}\\n")
    out.write(f"n_outliers={len(outliers)}\\n")
    for i, pv in enumerate(pve[:10]):
        out.write(f"pc{i+1}_variance_pct={pv:.2f}\\n")

print(f"Ancestry PCA: {len(outliers)} outliers flagged beyond ±{sd_thr} SD")
print(f"Covariate file: {n_study_written} study samples × {n_cov} PCs → pca_covariates.txt")
PYEOF

    # ── Scree plot: variance explained per PC ────────────────────────────────
    # Helps users decide how many PCs to include as GWAS/PRS covariates.
    # Shows individual and cumulative variance explained.
    if command -v Rscript &>/dev/null && [ -f pca.eigenval ]; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
evals <- scan("pca.eigenval", quiet=TRUE)
pve   <- evals / sum(evals) * 100
cpve  <- cumsum(pve)
df    <- data.frame(PC=seq_along(pve), PVE=pve, CPVE=cpve)
n_cov <- min(${params.n_pcs_covariates}, nrow(df))
p <- ggplot(df, aes(x=PC)) +
    geom_col(aes(y=PVE), fill="steelblue", alpha=0.8, width=0.7) +
    geom_line(aes(y=CPVE), colour="darkred", linewidth=0.8) +
    geom_point(aes(y=CPVE), colour="darkred", size=1.5) +
    geom_vline(xintercept=n_cov + 0.5, linetype="dashed", colour="grey40") +
    annotate("text", x=n_cov + 0.7, y=max(df$CPVE) * 0.5,
             label=paste0("n_pcs_covariates\n= ", n_cov),
             hjust=0, size=3, colour="grey40") +
    scale_x_continuous(breaks=seq_len(nrow(df))) +
    labs(title="PCA scree plot — variance explained per PC",
         subtitle="Bars: per-PC variance; red line: cumulative; dashed: covariate cutoff",
         x="Principal component", y="Variance explained (%)") +
    theme_classic() +
    theme(axis.text.x=element_text(size=7))
ggsave("pca_scree.png", p, width=9, height=5)
RSCRIPT
    fi

    # ── PCA scatter plot ──────────────────────────────────────────────────────
    # Priority: (1) population-labelled with HapMap info, (2) study vs reference,
    # (3) included vs outlier (no reference panel)
    if command -v Rscript &>/dev/null; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
pcs <- read.table("pca.eigenvec", header=FALSE)
colnames(pcs)[1:2] <- c("FID","IID")
colnames(pcs)[3:ncol(pcs)] <- paste0("PC", seq_len(ncol(pcs)-2))

study_fam <- read.table("${fam}", header=FALSE)
study_key <- paste(study_fam$V1, study_fam$V2)

hapmap_file <- "${hapmap_path}"
ref_merged  <- file.exists("merged_study_ref.fam")

if (ref_merged && nchar(hapmap_file) > 0 && file.exists(hapmap_file)) {
    # Best case: colour by HapMap population + mark study samples
    hapmap <- read.table(hapmap_file, header=TRUE)
    pcs <- merge(pcs, hapmap[, c("IID","population")], by="IID", all.x=TRUE)
    pcs$group <- ifelse(paste(pcs$FID, pcs$IID) %in% study_key, "Study",
                        ifelse(!is.na(pcs$population), pcs$population, "Reference"))
    pop_colours <- c(
        "Study"="steelblue",
        "CEU"="#E41A1C", "CHB"="#FF7F00", "JPT"="#984EA3",
        "YRI"="#4DAF4A", "ASW"="#A65628", "GIH"="#F781BF",
        "LWK"="#999999", "MEX"="#377EB8", "MKK"="#FFFF33",
        "TSI"="#00CED1", "CHD"="#8B4513"
    )
    colours  <- pop_colours[names(pop_colours) %in% unique(pcs$group)]
    subtitle <- "Coloured by HapMap reference population; blue = study samples"

} else if (ref_merged) {
    # Reference merged but no population info — study vs reference
    pcs$group <- ifelse(paste(pcs$FID, pcs$IID) %in% study_key, "Study", "Reference")
    colours   <- c("Study"="steelblue", "Reference"="grey60")
    subtitle  <- "Grey: reference panel samples"

} else {
    # No reference panel — colour ancestry outliers
    outliers <- tryCatch(
        read.table("ancestry_outliers.txt", header=FALSE, col.names=c("FID","IID")),
        error=function(e) data.frame(FID=character(0), IID=character(0))
    )
    pcs$group <- "Included"
    if (nrow(outliers) > 0) {
        key <- paste(outliers$FID, outliers$IID)
        pcs$group[paste(pcs$FID, pcs$IID) %in% key] <- "Outlier"
    }
    colours  <- c("Included"="steelblue", "Outlier"="red")
    subtitle <- paste0("Outlier threshold: ±${params.pca_outlier_sd} SD on any PC")
}

# Draw study samples on top so they're not hidden under reference points
pcs <- pcs[order(pcs$group != "Study"), ]

p <- ggplot(pcs, aes(x=PC1, y=PC2, colour=group)) +
    geom_point(alpha=0.6, size=ifelse(pcs$group == "Study", 1.2, 0.7)) +
    scale_colour_manual(values=colours) +
    labs(title="Ancestry PCA — PC1 vs PC2", subtitle=subtitle, colour="Group") +
    theme_classic() +
    guides(colour=guide_legend(override.aes=list(size=3, alpha=1)))
ggsave("pca_plot.png", p, width=9, height=6)
RSCRIPT
    fi
    """
}
