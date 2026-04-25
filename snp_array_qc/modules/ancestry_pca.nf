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

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val   reference_panel   // path prefix or empty list

    output:
    path "ancestry_outliers.txt",  emit: outlier_samples
    path "pca_summary.txt",        emit: summary
    path "pca.eigenvec",           emit: eigenvec
    path "pca.eigenval",           emit: eigenval

    script:
    def prefix = "${meta.id}"
    def has_ref = reference_panel instanceof List ? false : (reference_panel as String).length() > 0
    """
    # ── LD pruning ────────────────────────────────────────────────────────────
    # Use common, LD-independent SNPs to ensure PCA reflects ancestry, not LD
    plink \\
        --bfile ${bed.baseName} \\
        --maf 0.05 \\
        --geno 0.01 \\
        --hwe 0.001 \\
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

    # ── Identify ancestry outliers ────────────────────────────────────────────
    python3 - << 'PYEOF'
import math, sys

data = []
with open("pca.eigenvec") as fh:
    for line in fh:
        parts = line.split()
        fid, iid = parts[0], parts[1]
        pcs = list(map(float, parts[2:]))
        data.append((fid, iid, pcs))

n_pcs_computed = len(data[0][2]) if data else 0

# Compute mean and SD for each PC across all samples
means = [sum(d[2][i] for d in data) / len(data) for i in range(n_pcs_computed)]
sds   = [math.sqrt(sum((d[2][i] - means[i])**2 for d in data) / (len(data) - 1))
         for i in range(n_pcs_computed)]

sd_thr = ${params.pca_outlier_sd}
outliers = []
for fid, iid, pcs in data:
    for i, (pc, m, s) in enumerate(zip(pcs, means, sds)):
        if s > 0 and abs(pc - m) > sd_thr * s:
            outliers.append((fid, iid))
            break

with open("ancestry_outliers.txt", "w") as out:
    for fid, iid in outliers:
        out.write(f"{fid}\t{iid}\n")

# Read eigenvalues for variance explained
eigenvals = []
try:
    with open("pca.eigenval") as fh:
        eigenvals = [float(l.strip()) for l in fh if l.strip()]
    total_var = sum(eigenvals)
    pve = [v / total_var * 100 for v in eigenvals]
except Exception:
    pve = []

with open("pca_summary.txt", "w") as out:
    out.write(f"step=ancestry_pca\n")
    out.write(f"dataset=${meta.id}\n")
    out.write(f"pca_outlier_sd=${params.pca_outlier_sd}\n")
    out.write(f"n_pcs={n_pcs_computed}\n")
    out.write(f"n_samples={len(data)}\n")
    out.write(f"n_outliers={len(outliers)}\n")
    if pve:
        out.write(f"pc1_variance_pct={pve[0]:.2f}\n")
        out.write(f"pc2_variance_pct={pve[1]:.2f}\n")

print(f"Ancestry PCA: {len(outliers)} outliers flagged beyond ±{sd_thr} SD")
PYEOF

    # ── Optional: PCA scatter plot ────────────────────────────────────────────
    if command -v Rscript &>/dev/null; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
pcs <- read.table("pca.eigenvec", header=FALSE)
colnames(pcs)[1:2] <- c("FID","IID")
colnames(pcs)[3:ncol(pcs)] <- paste0("PC", seq_len(ncol(pcs)-2))
outliers <- tryCatch(read.table("ancestry_outliers.txt", header=FALSE,
                                col.names=c("FID","IID")), error=function(e) NULL)
pcs\$status <- "Included"
if (!is.null(outliers) && nrow(outliers) > 0) {
    key <- paste(outliers\$FID, outliers\$IID)
    pcs\$status[paste(pcs\$FID, pcs\$IID) %in% key] <- "Outlier"
}
p <- ggplot(pcs, aes(x=PC1, y=PC2, colour=status)) +
    geom_point(alpha=0.6, size=1) +
    scale_colour_manual(values=c("Included"="steelblue", "Outlier"="red")) +
    labs(title="Ancestry PCA — PC1 vs PC2",
         colour="Sample status") +
    theme_classic()
ggsave("pca_plot.png", p, width=8, height=6)
RSCRIPT
    fi
    """
}
