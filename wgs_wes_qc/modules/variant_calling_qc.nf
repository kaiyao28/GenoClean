/*
================================================================================
  MODULE: VARIANT_CALLING_QC
================================================================================
  Purpose:
    Compute variant-level and sample-level QC statistics from a VCF/gVCF file:
    SNP count, indel count, Ti/Tv ratio, het/hom ratio, singleton count,
    missingness, depth and GQ distributions.

  Why this step:
    Ti/Tv ratio is a classical QC metric: for WGS ~2.0–2.1 (whole genome),
    for WES ~2.8–3.0 (enriched for exonic variants). Values outside these
    ranges suggest variant calling errors or contamination.
    Singleton count is elevated in contaminated samples or those with alignment
    artefacts. Sample-level het/hom ratio can flag sample swaps.

  Output:
    - variant_stats.txt         : bcftools stats summary
    - sample_variant_stats.tsv  : per-sample SNP/indel/Ti/Tv/het/hom
    - variant_calling_qc_summary.txt : parsed summary for final report
================================================================================
*/

process VARIANT_CALLING_QC {
    label 'process_medium'
    publishDir "${params.outdir}/variant_calling_qc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)
    path reference_fasta

    output:
    path "variant_stats.txt",              emit: stats
    path "sample_variant_stats.tsv",       emit: sample_stats
    path "variant_calling_qc_summary.txt", emit: summary

    script:
    """
    # ── bcftools stats ────────────────────────────────────────────────────────
    bcftools stats \\
        --threads ${task.cpus} \\
        --fasta-ref ${reference_fasta} \\
        ${vcf} > variant_stats.txt

    # ── Per-sample statistics ─────────────────────────────────────────────────
    # bcftools stats -s - produces per-sample sections (PSC lines)
    bcftools stats \\
        --threads ${task.cpus} \\
        --fasta-ref ${reference_fasta} \\
        --samples - \\
        ${vcf} | grep "^PSC" > psc_lines.txt || true

    # ── Parse and report ──────────────────────────────────────────────────────
    python3 - << 'PYEOF'
import re

def parse_bcftools_stats(fname):
    metrics = {}
    try:
        with open(fname) as fh:
            for line in fh:
                if line.startswith("SN"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        key = parts[2].rstrip(":").replace(" ", "_")
                        metrics[key] = parts[3]
                elif line.startswith("TSTV"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 6:
                        metrics["ts_tv_ratio"] = parts[4]
    except Exception as e:
        print(f"Warning: {e}")
    return metrics

m = parse_bcftools_stats("variant_stats.txt")

n_snps       = m.get("number_of_SNPs", "NA")
n_indels     = m.get("number_of_indels", "NA")
ts_tv        = m.get("ts_tv_ratio", "NA")
n_multiallelic = m.get("number_of_multiallelic_sites", "NA")

# Parse PSC lines for per-sample het/hom
with open("sample_variant_stats.tsv", "w") as out:
    out.write("sample\tn_snps\tn_indels\tts_tv\tn_hom\tn_het\\n")
    try:
        with open("psc_lines.txt") as fh:
            for line in fh:
                parts = line.strip().split("\t")
                # PSC: id ts_tv n_hom n_het n_ts n_tv n_indel ...
                if len(parts) >= 14:
                    sample = parts[2]
                    n_hom  = parts[4]
                    n_het  = parts[5]
                    ts_tv_s = parts[6]
                    out.write(f"{sample}\tNA\tNA\t{ts_tv_s}\t{n_hom}\t{n_het}\\n")
    except Exception:
        out.write(f"${meta.id}\t{n_snps}\t{n_indels}\t{ts_tv}\tNA\tNA\\n")

with open("variant_calling_qc_summary.txt", "w") as out:
    out.write(f"step=variant_calling_qc\\n")
    out.write(f"dataset=${meta.id}\\n")
    out.write(f"n_snps={n_snps}\\n")
    out.write(f"n_indels={n_indels}\\n")
    out.write(f"ts_tv_ratio={ts_tv}\\n")
    out.write(f"n_multiallelic={n_multiallelic}\\n")

print(f"Variant calling QC: {n_snps} SNPs, {n_indels} indels, Ti/Tv={ts_tv}")
PYEOF
    """
}
