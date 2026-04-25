/*
================================================================================
  MODULE: ANCESTRY_PCA (WGS/WES)
================================================================================
  Purpose:
    Perform PCA and ancestry outlier detection on a pre-filtered VCF produced
    by PCA_VARIANT_SELECTION. This module does NOT select variants — that is
    the responsibility of the upstream PCA_VARIANT_SELECTION module.

  Design distinction from SNP array PCA:
    SNP ARRAY pipeline:
      - LD prune → PCA (all in one module)
      - Variant selection not needed: arrays genotype common variants by design

    WGS/WES pipeline:
      - PCA_VARIANT_SELECTION → LD prune → PCA (split across two modules)
      - Variant selection MUST happen first to exclude rare variants, indels,
        and high-missingness sites that distort PCA in sequencing data

  This module receives the clean common-SNP VCF and performs:
    1. LD pruning
    2. PLINK2 PCA
    3. Outlier detection (±params.pca_outlier_sd SD on any PC)
    4. Optional: projection / merge with reference panel
    5. PC scatter plot

  Shared parameters (same as SNP array PCA):
    params.n_pcs         = 20
    params.pca_outlier_sd = 6
    params.ld_window     = 200
    params.ld_step       = 50
    params.ld_r2         = 0.2
    params.reference_panel = null  (optional)

  How to disable:
    params.run_ancestry_pca_wgs = false

  Output:
    - pca.eigenvec               : sample PC scores
    - pca.eigenval               : variance explained per PC
    - ancestry_outliers.txt      : FID/IID of samples flagged as outliers
    - ancestry_pca_wgs_summary.txt : summary for the final report
    - ancestry_pca_wgs_plot.png  : PC1 vs PC2 scatter plot (if R available)
================================================================================
*/

process ANCESTRY_PCA {
    label 'process_high'
    publishDir "${params.outdir}/ancestry_pca", mode: params.publish_dir_mode

    input:
    // Input is a pre-filtered VCF from PCA_VARIANT_SELECTION (common biallelic SNPs only)
    tuple val(meta), path(vcf), path(tbi)
    val   reference_panel   // optional reference panel path prefix or []

    output:
    path "pca.eigenvec",                 emit: eigenvec
    path "pca.eigenval",                 emit: eigenval
    path "ancestry_outliers.txt",        emit: outliers
    path "ancestry_pca_wgs_summary.txt", emit: summary

    script:
    def has_ref = !(reference_panel instanceof List) && (reference_panel as String).length() > 0
    def ref_used = has_ref ? 'yes' : 'no'
    """
    n_variants=\$(bcftools view --no-header ${vcf} | wc -l)
    echo "ANCESTRY PCA: received \${n_variants} pre-filtered variants from PCA_VARIANT_SELECTION"

    # ── LD pruning ────────────────────────────────────────────────────────────
    # The input VCF is already filtered to common SNPs, but LD pruning is still
    # needed to ensure PCs reflect ancestry rather than local LD blocks.
    plink2 \\
        --vcf ${vcf} \\
        --indep-pairwise ${params.ld_window}kb ${params.ld_step} ${params.ld_r2} \\
        --out prune_pca \\
        --allow-extra-chr

    n_pruned=\$(wc -l < prune_pca.prune.in)
    echo "  After LD pruning: \${n_pruned} variants retained for PCA"

    # ── PCA ───────────────────────────────────────────────────────────────────
    ${has_ref ? """
    # Merge with reference panel for ancestry assignment
    # Extract the pruned variants from the reference panel
    plink2 \\
        --bfile ${reference_panel} \\
        --extract prune_pca.prune.in \\
        --export vcf bgz \\
        --out ref_subset \\
        --allow-extra-chr || true

    if [ -f ref_subset.vcf.gz ]; then
        bcftools index --tbi ref_subset.vcf.gz

        # Merge study VCF with reference subset
        bcftools merge \\
            --threads ${task.cpus} \\
            ${vcf} ref_subset.vcf.gz \\
            -O z -o merged_study_ref.vcf.gz || {
            echo "Merge failed — running PCA on study data only"
            plink2 \\
                --vcf ${vcf} \\
                --extract prune_pca.prune.in \\
                --pca ${params.n_pcs} \\
                --out pca \\
                --allow-extra-chr
        }

        if [ -f merged_study_ref.vcf.gz ]; then
            bcftools index --tbi merged_study_ref.vcf.gz
            plink2 \\
                --vcf merged_study_ref.vcf.gz \\
                --extract prune_pca.prune.in \\
                --pca ${params.n_pcs} \\
                --out pca \\
                --allow-extra-chr
        fi
    else
        echo "Reference panel conversion failed — running PCA on study data only"
        plink2 \\
            --vcf ${vcf} \\
            --extract prune_pca.prune.in \\
            --pca ${params.n_pcs} \\
            --out pca \\
            --allow-extra-chr
    fi
    """ : """
    # No reference panel — PCA on study data only.
    # PC interpretation will reflect internal structure only;
    # ancestry labels cannot be assigned without a reference.
    plink2 \\
        --vcf ${vcf} \\
        --extract prune_pca.prune.in \\
        --pca ${params.n_pcs} \\
        --out pca \\
        --allow-extra-chr
    """}

    # ── Outlier detection ─────────────────────────────────────────────────────
    python3 - << 'PYEOF'
import math

data = []
try:
    with open("pca.eigenvec") as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.split()
            fid, iid = parts[0], parts[1]
            pcs = list(map(float, parts[2:]))
            data.append((fid, iid, pcs))
except Exception as e:
    print(f"Warning: {e}")

if not data:
    with open("ancestry_outliers.txt", "w") as out:
        pass
    with open("ancestry_pca_wgs_summary.txt", "w") as out:
        out.write("step=ancestry_pca_wgs\ndataset=${meta.id}\nn_samples=0\nn_outliers=0\\n")
    exit()

n_pcs_computed = len(data[0][2])
means = [sum(d[2][i] for d in data) / len(data) for i in range(n_pcs_computed)]
sds   = [math.sqrt(sum((d[2][i] - means[i])**2 for d in data) / max(len(data) - 1, 1))
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
        out.write(f"{fid}\t{iid}\\n")

eigenvals = []
pve = []
try:
    with open("pca.eigenval") as fh:
        eigenvals = [float(l.strip()) for l in fh if l.strip()]
    total_var = sum(eigenvals)
    pve = [v / total_var * 100 for v in eigenvals] if total_var > 0 else []
except Exception:
    pass

with open("ancestry_pca_wgs_summary.txt", "w") as out:
    out.write(f"step=ancestry_pca_wgs\\n")
    out.write(f"dataset=${meta.id}\\n")
    out.write(f"mode=${meta.mode}\\n")
    out.write(f"pca_outlier_sd=${params.pca_outlier_sd}\\n")
    out.write(f"n_pcs_computed={n_pcs_computed}\\n")
    out.write(f"n_samples={len(data)}\\n")
    out.write(f"n_outliers={len(outliers)}\\n")
    if pve:
        out.write(f"pc1_variance_pct={pve[0]:.2f}\\n")
        out.write(f"pc2_variance_pct={pve[1]:.2f}\\n" if len(pve) > 1 else "")
    out.write("reference_panel_used=${ref_used}\\n")

print(f"Ancestry PCA WGS/WES: {len(outliers)} outliers flagged beyond ±{sd_thr} SD "
      f"({len(data)} samples, {n_pcs_computed} PCs)")
PYEOF

    # ── PC scatter plot ───────────────────────────────────────────────────────
    if command -v Rscript &>/dev/null && [ -f pca.eigenvec ]; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
pcs <- read.table("pca.eigenvec", header=TRUE)
if (ncol(pcs) < 4) quit(status=0)
colnames(pcs)[1:2] <- c("FID","IID")
pcs\$status <- "Included"
outliers <- tryCatch(
    read.table("ancestry_outliers.txt", col.names=c("FID","IID")),
    error=function(e) NULL
)
if (!is.null(outliers) && nrow(outliers) > 0) {
    pcs\$status[paste(pcs\$FID, pcs\$IID) %in% paste(outliers\$FID, outliers\$IID)] <- "Outlier"
}
eigenvals <- tryCatch(scan("pca.eigenval", quiet=TRUE), error=function(e) NULL)
xlab <- if (!is.null(eigenvals) && length(eigenvals) >= 1)
    sprintf("PC1 (%.1f%%)", eigenvals[1]/sum(eigenvals)*100) else "PC1"
ylab <- if (!is.null(eigenvals) && length(eigenvals) >= 2)
    sprintf("PC2 (%.1f%%)", eigenvals[2]/sum(eigenvals)*100) else "PC2"
p <- ggplot(pcs, aes(x=pcs[[3]], y=pcs[[4]], colour=status)) +
    geom_point(alpha=0.6, size=1.5) +
    scale_colour_manual(values=c("Included"="steelblue", "Outlier"="red")) +
    labs(title=paste0("Ancestry PCA — WGS/WES — ${meta.id}"),
         subtitle=sprintf("±${params.pca_outlier_sd} SD outlier threshold | %d samples",
                          nrow(pcs)),
         x=xlab, y=ylab, colour="Status") +
    theme_classic()
ggsave("ancestry_pca_wgs_plot.png", p, width=8, height=6, dpi=150)
RSCRIPT
    fi
    """
}
