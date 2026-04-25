/*
================================================================================
  MODULE: RELATEDNESS
================================================================================
  Purpose:
    Identify pairs of samples with unexpectedly high genetic relatedness (IBD).
    Cryptic relatedness violates the independence assumption of GWAS and inflates
    test statistics. One sample from each related pair is recommended for removal.

  Why this step:
    PLINK estimates PI_HAT (proportion of alleles shared IBD) for all pairs.
    Expected PI_HAT values:
      Monozygotic twins / duplicate samples : ~1.0
      First-degree relatives (parent-child, full siblings) : ~0.5
      Second-degree relatives (half-siblings, grandparent-grandchild) : ~0.25
      Third-degree relatives : ~0.125

    The default threshold 0.1875 sits between 2nd and 3rd degree, capturing
    all relatives up to and including 2nd degree. For population isolates or
    small cohorts, a more permissive threshold may be appropriate.

  Default threshold: params.relatedness_pi_hat = 0.1875

  How to change:
    nextflow run main.nf --relatedness_pi_hat 0.125   # include 3rd-degree

  How to disable:
    nextflow run main.nf --run_relatedness false

  Note:
    This module removes one member of each related pair using a greedy algorithm
    (the sample with more relatives is removed first). Both members of duplicate
    pairs (PI_HAT > 0.9) are flagged; the researcher should decide which to keep.

  Output:
    - relatedness.genome       : PLINK IBD output (related pairs only)
    - relatedness_remove.txt   : FID IID recommended for removal
    - relatedness_summary.txt  : counts
================================================================================
*/

process RELATEDNESS {
    label 'process_high'
    publishDir "${params.outdir}/qc_tables", mode: params.publish_dir_mode, pattern: "*.{genome,txt}"
    publishDir "${params.outdir}/qc_plots",  mode: params.publish_dir_mode, pattern: "*.png"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    path "relatedness_remove.txt", emit: related_samples
    path "relatedness_summary.txt", emit: summary
    path "relatedness.genome",      emit: genome

    script:
    def prefix = "${meta.id}"
    """
    # ── LD pruning: IBD estimation is more accurate on independent SNPs ────────
    plink \\
        --bfile ${bed.baseName} \\
        --maf 0.05 \\
        --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \\
        --out prune_rel \\
        --allow-no-sex

    # ── IBD estimation ────────────────────────────────────────────────────────
    # --genome outputs pairwise IBD estimates
    # --min filters to pairs exceeding the PI_HAT threshold (saves disk space)
    plink \\
        --bfile ${bed.baseName} \\
        --extract prune_rel.prune.in \\
        --genome \\
        --min ${params.relatedness_pi_hat} \\
        --out relatedness \\
        --allow-no-sex

    # ── Select one sample from each related pair for removal ─────────────────
    python3 - << 'PYEOF'
import sys
from collections import defaultdict

pairs = []
with open("relatedness.genome") as fh:
    next(fh)  # skip header
    for line in fh:
        parts = line.split()
        if len(parts) < 10:
            continue
        fid1, iid1 = parts[0], parts[1]
        fid2, iid2 = parts[2], parts[3]
        pi_hat = float(parts[9])
        pairs.append(((fid1, iid1), (fid2, iid2), pi_hat))

n_pairs = len(pairs)
n_duplicates = sum(1 for _, _, pi in pairs if pi > 0.9)

# Count how many related pairs each individual appears in
degree = defaultdict(int)
for s1, s2, pi in pairs:
    degree[s1] += 1
    degree[s2] += 1

# Greedy removal: iteratively remove the sample with the most connections
to_remove = set()
remaining_pairs = list(pairs)
while remaining_pairs:
    # Recount connections in surviving pairs
    count = defaultdict(int)
    for s1, s2, _ in remaining_pairs:
        if s1 not in to_remove and s2 not in to_remove:
            count[s1] += 1
            count[s2] += 1
    if not count:
        break
    # Remove the individual with the most connections
    victim = max(count, key=count.get)
    to_remove.add(victim)
    remaining_pairs = [(s1, s2, pi) for s1, s2, pi in remaining_pairs
                       if victim not in (s1, s2)]

with open("relatedness_remove.txt", "w") as out:
    for fid, iid in to_remove:
        out.write(f"{fid}\t{iid}\n")

with open("relatedness_summary.txt", "w") as out:
    out.write(f"step=relatedness\n")
    out.write(f"dataset=${meta.id}\n")
    out.write(f"pi_hat_threshold=${params.relatedness_pi_hat}\n")
    out.write(f"n_related_pairs={n_pairs}\n")
    out.write(f"n_duplicate_pairs={n_duplicates}\n")
    out.write(f"n_samples_removed={len(to_remove)}\n")

print(f"Relatedness: {n_pairs} pairs above PI_HAT {${params.relatedness_pi_hat}}, "
      f"{len(to_remove)} samples recommended for removal")
PYEOF

    # ── Optional: PI_HAT distribution plot ────────────────────────────────────
    if command -v Rscript &>/dev/null && [ -s relatedness.genome ]; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
df <- read.table("relatedness.genome", header=TRUE)
p <- ggplot(df, aes(x=PI_HAT)) +
    geom_histogram(bins=100, fill="steelblue", alpha=0.8) +
    geom_vline(xintercept=${params.relatedness_pi_hat}, linetype="dashed", colour="red") +
    labs(title="Pairwise IBD (PI_HAT) distribution",
         subtitle=paste0("Threshold: ", ${params.relatedness_pi_hat}),
         x="PI_HAT", y="Number of pairs") +
    theme_classic()
ggsave("relatedness_pi_hat.png", p, width=8, height=5)
RSCRIPT
    fi
    """
}
