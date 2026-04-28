/*
================================================================================
  MODULE: IMPUTATION_FILTER
================================================================================
  Purpose:
    Filter poorly imputed variants using an imputation quality score.
    Only relevant for post-imputation QC — leave disabled (default) when
    working with directly genotyped SNP array data.

  Supported formats (auto-detected from file header):
    Michigan / Minimac4  : tab-separated .info or .info.gz
                           Quality column : Rsq
                           SNP ID column  : SNP (format chr:pos:ref:alt)
    IMPUTE2 / IMPUTE5    : space-separated .info or .info.gz
                           Quality column : info (column 5)
                           SNP ID column  : rs_id (column 2)
    BEAGLE               : VCF (.vcf or .vcf.gz) with DR2 in INFO field
                           Quality column : DR2 (parsed from INFO)
                           SNP ID column  : VCF ID field (column 3)

  Default threshold: params.imputation_r2 = 0.3
    - 0.3 : Michigan's "acceptable quality" cutoff for common variants
    - 0.8 : stringent; recommended for rare variant analysis

  How to enable:
    nextflow run main.nf --run_imputation_filter true \
                         --info_file results/imputed/chr_all.info.gz

  How to change threshold:
    nextflow run main.nf --run_imputation_filter true \
                         --info_file data/imputed.info \
                         --imputation_r2 0.8

  Note on SNP ID matching:
    The IDs in the info file must match column 2 of the .bim file.
    Michigan produces chr:pos:ref:alt IDs — ensure your PLINK files use
    the same convention (e.g. converted with plink2 --set-all-var-ids).
    A low-match warning is printed if fewer than 10% of info IDs are found
    in the .bim file.

  Output:
    - Filtered PLINK files
    - imputation_fail.txt           : variant IDs below R2 threshold
    - imputation_filter_summary.txt : format detected, threshold, counts
    - imputation_r2_plot.png        : R2 score distribution (if R available)
================================================================================
*/

process IMPUTATION_FILTER {
    label 'process_medium'
    publishDir "${params.outdir}/qc_tables",    mode: params.publish_dir_mode, pattern: "*.txt"
    publishDir "${params.outdir}/qc_plots",     mode: params.publish_dir_mode, pattern: "*.png"
    publishDir "${params.outdir}/cleaned_data", mode: params.publish_dir_mode, pattern: "*.{bed,bim,fam}", enabled: params.keep_intermediate

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path info_file   // Michigan .info[.gz], IMPUTE2 .info[.gz], or BEAGLE VCF; [] if not provided

    output:
    tuple val(meta), path("${meta.id}_impfilt.bed"),
                     path("${meta.id}_impfilt.bim"),
                     path("${meta.id}_impfilt.fam"), emit: plink
    path "imputation_fail.txt",             emit: removed_variants
    path "imputation_filter_summary.txt",   emit: summary
    path "*.png",                           optional: true, emit: plots

    script:
    def prefix        = "${meta.id}"
    def has_info      = !(info_file instanceof List)
    def info_path_str = has_info ? "\"${info_file}\"" : "\"\""
    """
    n_var_before=\$(wc -l < ${bim})

    # ── Parse info file and build exclusion list ──────────────────────────────
    python3 - << 'PYEOF'
import gzip, sys, os

info_path = ${info_path_str}   # empty string when no info file provided
r2_thresh = float("${params.imputation_r2}")

def open_any(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

if not info_path:
    # No info file — write empty outputs and exit cleanly
    open("imputation_fail.txt",      "w").close()
    open("imputation_r2_scores.txt", "w").close()
    with open("imputation_filter_summary.txt", "w") as fh:
        fh.write("step=imputation_filter\\n")
        fh.write("dataset=${meta.id}\\n")
        fh.write("format=none\\n")
        fh.write("r2_threshold=${params.imputation_r2}\\n")
        fh.write("note=no_info_file_provided\\n")
    print("No info file provided — imputation filter pass-through")
    sys.exit(0)

# ── Read .bim IDs for overlap check ──────────────────────────────────────────
bim_ids = set()
with open("${bim}") as fh:
    for line in fh:
        p = line.split()
        if len(p) >= 2:
            bim_ids.add(p[1])

# ── Detect format from first non-comment line ─────────────────────────────────
with open_any(info_path) as fh:
    first_line = fh.readline()

header = first_line.strip()

if "Rsq" in header and "ALT_Frq" in header:
    fmt = "michigan"
elif "exp_freq_a1" in header:
    fmt = "impute2"
elif header.startswith("##") or header.startswith("#CHROM"):
    fmt = "beagle_vcf"
else:
    sys.exit(f"ERROR: Cannot detect info file format.\\nHeader: {header[:300]}")

print(f"Detected format: {fmt}")

# ── Parse scores ──────────────────────────────────────────────────────────────
fail_ids  = []
r2_scores = []   # list of (snp_id, score)
n_total   = 0

if fmt == "michigan":
    with open_any(info_path) as fh:
        cols = fh.readline().strip().split("\\t")
        try:
            snp_idx = cols.index("SNP")
            rsq_idx = cols.index("Rsq")
        except ValueError:
            sys.exit("ERROR: Michigan .info file missing SNP or Rsq column")
        for line in fh:
            p = line.strip().split("\\t")
            if len(p) <= max(snp_idx, rsq_idx):
                continue
            sid = p[snp_idx]
            try:
                score = float(p[rsq_idx])
            except ValueError:
                continue
            n_total += 1
            r2_scores.append((sid, score))
            if score < r2_thresh:
                fail_ids.append(sid)

elif fmt == "impute2":
    # Columns: snp_id  rs_id  position  exp_freq_a1  info  certainty  type ...
    with open_any(info_path) as fh:
        fh.readline()  # skip header
        for line in fh:
            p = line.strip().split()
            if len(p) < 5:
                continue
            sid = p[1]        # rs_id
            try:
                score = float(p[4])   # info score
            except ValueError:
                continue
            n_total += 1
            r2_scores.append((sid, score))
            if score < r2_thresh:
                fail_ids.append(sid)

elif fmt == "beagle_vcf":
    # DR2 embedded in VCF INFO field as DR2=<value>
    with open_any(info_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.strip().split("\\t")
            if len(p) < 8:
                continue
            sid        = p[2]    # VCF ID column
            info_field = p[7]
            for token in info_field.split(";"):
                if token.startswith("DR2="):
                    try:
                        score = float(token[4:])
                        n_total += 1
                        r2_scores.append((sid, score))
                        if score < r2_thresh:
                            fail_ids.append(sid)
                    except ValueError:
                        pass
                    break

# ── Overlap check ─────────────────────────────────────────────────────────────
info_snp_ids  = {sid for sid, _ in r2_scores}
n_bim_matched = len(info_snp_ids & bim_ids)
overlap_pct   = n_bim_matched / max(len(info_snp_ids), 1) * 100
if overlap_pct < 10:
    print(
        f"WARNING: Only {overlap_pct:.1f}% of info file IDs match .bim SNP IDs.\\n"
        f"Check ID format consistency — Michigan uses chr:pos:ref:alt, "
        f"IMPUTE2 uses rsIDs. Align with plink2 --set-all-var-ids if needed."
    )

# ── Write outputs ─────────────────────────────────────────────────────────────
with open("imputation_fail.txt", "w") as fh:
    for sid in fail_ids:
        fh.write(sid + "\\n")

with open("imputation_r2_scores.txt", "w") as fh:
    fh.write("SNP\\tR2\\n")
    for sid, score in r2_scores:
        fh.write(f"{sid}\\t{score}\\n")

n_fail = len(fail_ids)
with open("imputation_filter_summary.txt", "w") as fh:
    fh.write(f"step=imputation_filter\\n")
    fh.write(f"dataset=${meta.id}\\n")
    fh.write(f"format={fmt}\\n")
    fh.write(f"r2_threshold=${params.imputation_r2}\\n")
    fh.write(f"n_variants_in_info={n_total}\\n")
    fh.write(f"n_variants_failing={n_fail}\\n")
    fh.write(f"n_bim_ids_matched={n_bim_matched}\\n")
    fh.write(f"bim_id_overlap_pct={overlap_pct:.1f}\\n")

print(f"Imputation filter ({fmt}): {n_fail}/{n_total} variants below R2 {r2_thresh}")
print(f"Bim ID overlap: {n_bim_matched}/{len(info_snp_ids)} ({overlap_pct:.1f}%)")
PYEOF

    # ── Apply exclusion with PLINK ────────────────────────────────────────────
    n_fail_actual=\$(wc -l < imputation_fail.txt)

    if [ "\${n_fail_actual}" -gt 0 ]; then
        plink \\
            --bfile ${bed.baseName} \\
            --exclude imputation_fail.txt \\
            --make-bed \\
            --out ${prefix}_impfilt \\
            --allow-no-sex
    else
        plink \\
            --bfile ${bed.baseName} \\
            --make-bed \\
            --out ${prefix}_impfilt \\
            --allow-no-sex
    fi

    n_var_after=\$(wc -l < ${prefix}_impfilt.bim)
    n_removed=\$(( n_var_before - n_var_after ))

    echo "Imputation filter: removed \${n_removed} of \${n_var_before} variants (R2 < ${params.imputation_r2})"

    # ── R2 distribution plot ──────────────────────────────────────────────────
    if command -v Rscript &>/dev/null && [ -s imputation_r2_scores.txt ]; then
        Rscript - << 'RSCRIPT'
library(ggplot2)
df <- read.table("imputation_r2_scores.txt", header=TRUE, sep="\t")
df$R2 <- suppressWarnings(as.numeric(df$R2))
df <- df[!is.na(df$R2), ]
if (nrow(df) == 0) quit(save="no")
p <- ggplot(df, aes(x=R2)) +
    geom_histogram(bins=100, fill="steelblue", alpha=0.8) +
    geom_vline(xintercept=${params.imputation_r2}, linetype="dashed", colour="red") +
    scale_x_continuous(limits=c(0, 1)) +
    labs(title="Imputation quality (R²) distribution",
         subtitle=paste0("Format: auto-detected | Red line: threshold ${params.imputation_r2}"),
         x="Imputation R² / INFO score", y="Number of variants") +
    theme_classic()
ggsave("imputation_r2_plot.png", p, width=8, height=5)
RSCRIPT
    fi
    """
}
