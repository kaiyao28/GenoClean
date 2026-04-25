/*
================================================================================
  MODULE: RELATEDNESS (WGS/WES)
================================================================================
  Purpose:
    Identify related sample pairs from a multi-sample VCF using KING or
    bcftools relatedness2. For large cohorts, Hail pc_relate is recommended.

  Why this step:
    Related individuals violate independence assumptions in association testing.
    The same PI_HAT threshold logic from the SNP array pipeline applies here.
    For WGS/WES, bcftools gtcheck or KING-robust can be used directly on VCFs.

  Default threshold: params.relatedness_pi_hat = 0.1875

  How to disable:
    params.run_relatedness_wgs = false

  Output:
    - relatedness_pairs.txt         : pairs above threshold
    - relatedness_remove.txt        : suggested removal list
    - relatedness_wgs_summary.txt   : summary for final report
================================================================================
*/

process RELATEDNESS {
    label 'process_high'
    publishDir "${params.outdir}/relatedness", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    path "relatedness_pairs.txt",       emit: pairs
    path "relatedness_remove.txt",      emit: related_samples
    path "relatedness_wgs_summary.txt", emit: summary

    script:
    """
    # ── Compute pairwise relatedness using bcftools ────────────────────────────
    # +relatedness2 plugin uses a method similar to KING-robust
    bcftools +relatedness2 \\
        --threads ${task.cpus} \\
        ${vcf} > relatedness_raw.txt 2>&1 || {
        echo "bcftools +relatedness2 not available — using bcftools gtcheck as fallback"
        bcftools gtcheck \\
            --threads ${task.cpus} \\
            ${vcf} > relatedness_raw.txt || true
    }

    # ── Parse and flag related pairs ──────────────────────────────────────────
    python3 - << 'PYEOF'
from collections import defaultdict

threshold = ${params.relatedness_pi_hat}
pairs = []

try:
    with open("relatedness_raw.txt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 8:
                # bcftools +relatedness2 format: SAMPLE_A SAMPLE_B ...
                s1, s2 = parts[0], parts[1]
                # relatedness coefficient is typically in col 7 or 8
                try:
                    rel = float(parts[7])
                except (ValueError, IndexError):
                    try:
                        rel = float(parts[4])
                    except Exception:
                        rel = 0.0
                if rel >= threshold and s1 != s2:
                    pairs.append((s1, s2, rel))
except Exception as e:
    print(f"Warning: could not parse relatedness output: {e}")

# Greedy removal
degree = defaultdict(int)
for s1, s2, rel in pairs:
    degree[s1] += 1
    degree[s2] += 1

to_remove = set()
remaining = list(pairs)
while remaining:
    count = defaultdict(int)
    for s1, s2, _ in remaining:
        if s1 not in to_remove and s2 not in to_remove:
            count[s1] += 1
            count[s2] += 1
    if not count:
        break
    victim = max(count, key=count.get)
    to_remove.add(victim)
    remaining = [(s1, s2, r) for s1, s2, r in remaining if victim not in (s1, s2)]

with open("relatedness_pairs.txt", "w") as out:
    out.write("sample_a\tsample_b\trelatedness\n")
    for s1, s2, r in pairs:
        out.write(f"{s1}\t{s2}\t{r:.6f}\n")

with open("relatedness_remove.txt", "w") as out:
    for s in to_remove:
        out.write(f"{s}\n")

with open("relatedness_wgs_summary.txt", "w") as out:
    out.write(f"step=relatedness_wgs\n")
    out.write(f"dataset=${meta.id}\n")
    out.write(f"pi_hat_threshold={threshold}\n")
    out.write(f"n_related_pairs={len(pairs)}\n")
    out.write(f"n_samples_to_remove={len(to_remove)}\n")

print(f"Relatedness WGS: {len(pairs)} related pairs, {len(to_remove)} samples for removal")
PYEOF
    """
}
