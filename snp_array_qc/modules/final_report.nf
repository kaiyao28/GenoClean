/*
================================================================================
  MODULE: FINAL_REPORT (SNP Array)
================================================================================
  Purpose:
    Aggregate per-step SNP-array QC summaries into one HTML report and TSVs.
================================================================================
*/

process FINAL_REPORT {
    label 'process_report'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path summary_files
    path excluded_samples
    path excluded_variants
    path plot_files      // collected PNGs from all QC modules; may be empty list
    path qc_data_files   // raw per-sample QC data files (imiss/het/sexcheck/removal lists)
    val qc_scope

    output:
    path "final_report.html",       emit: report
    path "qc_report.pdf",           optional: true, emit: pdf
    path "qc_attrition_table.tsv",  emit: attrition
    path "qc_thresholds.tsv",       emit: thresholds
    path "qc_per_sample.tsv",       optional: true, emit: per_sample
    path "qc_per_batch.tsv",        optional: true, emit: per_batch

    script:
    """
    python3 - << 'PYEOF'
import base64
import datetime
import glob
import html
import math
import os
from collections import defaultdict

summaries = []
for fname in sorted(glob.glob("*summary*.txt")):
    data = {"file": fname}
    try:
        with open(fname) as fh:
            for line in fh:
                line = line.strip()
                if "=" in line:
                    k, v = line.split("=", 1)
                    data[k.strip()] = v.strip()
    except Exception as exc:
        data["parse_error"] = str(exc)
    summaries.append(data)

def count_lines(path):
    try:
        with open(path) as fh:
            return sum(1 for line in fh if line.strip())
    except Exception:
        return 0

final_samples = count_lines("${fam}")
final_variants = count_lines("${bim}")
excluded_sample_count = count_lines("${excluded_samples}")
excluded_variant_count = count_lines("${excluded_variants}")

variant_steps = {"duplicate_variant_check", "variant_callrate", "hwe_filter", "maf_filter"}
sample_steps = {"sample_callrate", "sex_check", "heterozygosity", "relatedness", "ancestry_pca"}

with open("qc_attrition_table.tsv", "w") as out:
    out.write("phase\\tstep\\tmetric\\tvalue\\n")
    for data in summaries:
        step = data.get("step", data["file"].replace(".txt", ""))
        phase = "variant_level_qc" if step in variant_steps else "sample_level_qc" if step in sample_steps else "input_check"
        for key, value in sorted(data.items()):
            if key in ("file", "step"):
                continue
            out.write(f"{phase}\\t{step}\\t{key}\\t{value}\\n")
    out.write(f"final\\tcleaned_data\\tn_samples\\t{final_samples}\\n")
    out.write(f"final\\tcleaned_data\\tn_variants\\t{final_variants}\\n")
    out.write(f"final\\texclusions\\tn_excluded_samples\\t{excluded_sample_count}\\n")
    out.write(f"final\\texclusions\\tn_excluded_variants\\t{excluded_variant_count}\\n")

thresholds = [
    ("chroms", "${params.chroms}"),
    ("sample_qc_scope", "${qc_scope}"),
    ("run_variant_qc", "${params.run_variant_qc}"),
    ("run_sample_qc", "${params.run_sample_qc}"),
    ("sample_missingness", "${params.sample_missingness}"),
    ("variant_missingness", "${params.variant_missingness}"),
    ("cc_miss_p", "${params.cc_miss_p}"),
    ("hwe_p_autosomes", "${params.hwe_p}"),
    ("hwe_p_chrx", "${params.hwe_p_chrx}"),
    ("maf", "${params.maf}"),
    ("heterozygosity_sd", "${params.heterozygosity_sd}"),
    ("relatedness_pi_hat", "${params.relatedness_pi_hat}"),
    ("sex_check_f_lower_female", "${params.sex_check_f_lower_female}"),
    ("sex_check_f_upper_male", "${params.sex_check_f_upper_male}"),
    ("pca_outlier_sd", "${params.pca_outlier_sd}"),
    ("ld_regions", "${params.ld_regions ?: 'not provided'}"),
    ("reference_panel", "${params.reference_panel ?: 'not provided'}"),
]

with open("qc_thresholds.tsv", "w") as out:
    out.write("parameter\\tvalue\\n")
    for key, value in thresholds:
        out.write(f"{key}\\t{value}\\n")

# ── Per-sample QC table ───────────────────────────────────────────────────────
def read_id_set(path):
    ids = set()
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 2:
                    ids.add((parts[0], parts[1]))
    except Exception:
        pass
    return ids

def extract_batch(fid):
    parts = fid.split("-")
    return "-".join(parts[:2]) if len(parts) >= 3 else fid

imiss_data = {}
if os.path.exists("sample_missingness.imiss"):
    with open("sample_missingness.imiss") as fh:
        next(fh)
        for line in fh:
            p = line.split()
            if len(p) >= 6:
                imiss_data[(p[0], p[1])] = float(p[5])

het_data = {}
if os.path.exists("heterozygosity.het"):
    with open("heterozygosity.het") as fh:
        next(fh)
        for line in fh:
            p = line.split()
            if len(p) >= 5:
                o_hom, n_nm = float(p[2]), float(p[4])
                het_data[(p[0], p[1])] = (n_nm - o_hom) / n_nm if n_nm > 0 else float("nan")

sex_data = {}
if os.path.exists("sex_check.sexcheck"):
    with open("sex_check.sexcheck") as fh:
        next(fh)
        for line in fh:
            p = line.split()
            if len(p) >= 6:
                sex_data[(p[0], p[1])] = (p[2], p[3], p[4], p[5])

smiss_removed  = read_id_set("sample_callrate_removed.txt")
het_outliers   = read_id_set("heterozygosity_outliers.txt")
sex_discordant = read_id_set("sex_discordant.txt")
related_remove = read_id_set("relatedness_remove.txt")
anc_outliers   = read_id_set("ancestry_outliers.txt")

seen, all_samples = set(), []
for k in list(imiss_data) + list(het_data) + list(sex_data):
    if k not in seen:
        seen.add(k); all_samples.append(k)

ps_header = ["FID","IID","F_MISS","HET_RATE","X_F_STAT","PEDSEX","SNPSEX","SEX_STATUS",
             "SMISS_REMOVE","HET_OUTLIER","SEX_DISCORDANT","RELATED_REMOVE","ANCESTRY_OUTLIER",
             "QC_FAIL","BATCH"]
sample_rows_data = []
for (fid, iid) in all_samples:
    f_miss   = imiss_data.get((fid, iid), float("nan"))
    het_rate = het_data.get((fid, iid), float("nan"))
    sx       = sex_data.get((fid, iid), ("NA","NA","NA","NA"))
    sm_rm  = "YES" if (fid,iid) in smiss_removed   else "NO"
    ht_out = "YES" if (fid,iid) in het_outliers    else "NO"
    sx_dis = "YES" if (fid,iid) in sex_discordant  else "NO"
    rel_rm = "YES" if (fid,iid) in related_remove  else "NO"
    anc_out= "YES" if (fid,iid) in anc_outliers    else "NO"
    qc_fail= "YES" if any(x=="YES" for x in [sm_rm,ht_out,sx_dis,rel_rm,anc_out]) else "NO"
    batch  = extract_batch(fid)
    row = [fid, iid,
           f"{f_miss:.5f}" if not math.isnan(f_miss) else "NA",
           f"{het_rate:.5f}" if not math.isnan(het_rate) else "NA",
           sx[3], sx[0], sx[1], sx[2],
           sm_rm, ht_out, sx_dis, rel_rm, anc_out, qc_fail, batch]
    sample_rows_data.append(row)

if sample_rows_data:
    with open("qc_per_sample.tsv", "w") as out:
        out.write("\\t".join(ps_header) + "\\n")
        for row in sample_rows_data:
            out.write("\\t".join(str(x) for x in row) + "\\n")

batch_n    = defaultdict(int)
batch_fail = defaultdict(int)
for row in sample_rows_data:
    b = row[14]
    batch_n[b] += 1
    if row[13] == "YES":
        batch_fail[b] += 1

batch_rows_data = []
for b in sorted(batch_n):
    n = batch_n[b]; r = batch_fail[b]
    pct = f"{(n-r)/n*100:.1f}" if n > 0 else "NA"
    batch_rows_data.append([b, n, r, n-r, pct])

if batch_rows_data:
    with open("qc_per_batch.tsv", "w") as out:
        out.write("BATCH\\tN_SAMPLES\\tN_REMOVED\\tN_PASS\\tPASS_RATE_PCT\\n")
        for row in batch_rows_data:
            out.write("\\t".join(str(x) for x in row) + "\\n")

print(f"Per-sample table: {len(sample_rows_data)} samples")
print(f"Per-batch table: {len(batch_rows_data)} batches")

summary_rows = []
with open("qc_attrition_table.tsv") as fh:
    next(fh)
    for line in fh:
        phase, step, metric, value = line.rstrip("\\n").split("\\t", 3)
        summary_rows.append(
            f"<tr><td>{html.escape(phase)}</td><td>{html.escape(step)}</td>"
            f"<td>{html.escape(metric)}</td><td>{html.escape(value)}</td></tr>"
        )

threshold_rows = [
    f"<tr><td>{html.escape(k)}</td><td>{html.escape(v)}</td></tr>"
    for k, v in thresholds
]

# ── Per-batch HTML table ──────────────────────────────────────────────────────
batch_html = ""
if batch_rows_data:
    brows = "".join(
        f"<tr><td>{html.escape(str(r[0]))}</td><td>{r[1]}</td><td>{r[2]}</td>"
        f"<td>{r[3]}</td><td>{r[4]}%</td></tr>"
        for r in batch_rows_data
    )
    batch_html = f"""<h2>Per-batch QC Summary</h2>
<table>
<tr><th>Batch</th><th>Samples</th><th>Removed</th><th>Pass</th><th>Pass %</th></tr>
{brows}
</table>"""

# ── Per-sample HTML table ─────────────────────────────────────────────────────
per_sample_html = ""
if sample_rows_data:
    def flag_td(val):
        cls = "yes" if val == "YES" else "no"
        return f'<td class="{cls}">{val}</td>'
    srows = "".join(
        f"<tr><td>{html.escape(r[0])}</td><td>{html.escape(r[1])}</td>"
        f"<td>{r[2]}</td><td>{r[3]}</td><td>{r[4]}</td>"
        f"<td>{r[5]}</td><td>{r[6]}</td><td>{r[7]}</td>"
        + flag_td(r[8]) + flag_td(r[9]) + flag_td(r[10])
        + flag_td(r[11]) + flag_td(r[12]) + flag_td(r[13])
        + f"<td>{html.escape(r[14])}</td></tr>"
        for r in sample_rows_data
    )
    per_sample_html = f"""<h2>Per-sample QC Table</h2>
<p>Full table available as <code>qc_per_sample.tsv</code>. Flags in red indicate QC failure at that step.</p>
<div class="scroll-table">
<table>
<tr><th>FID</th><th>IID</th><th>F_MISS</th><th>HET_RATE</th><th>X_F_STAT</th>
<th>PEDSEX</th><th>SNPSEX</th><th>SEX_STATUS</th>
<th>SMISS</th><th>HET_OUT</th><th>SEX_DISC</th><th>RELATED</th><th>ANC_OUT</th><th>QC_FAIL</th><th>BATCH</th></tr>
{srows}
</table>
</div>"""

# ── Collect and base64-encode all QC plots ────────────────────────────────
plot_order = [
    "miss_het_scatter",
    "heterozygosity_plot",
    "sex_check_F_stat",
    "vmiss_plot",
    "cc_miss_plot",
    "maf_plot",
    "hwe_plot",
    "imputation_r2_plot",
    "relatedness_pi_hat",
    "pca_scree",
    "pca_plot",
]

def plot_sort_key(path):
    name = os.path.splitext(os.path.basename(path))[0]
    try:
        return plot_order.index(name)
    except ValueError:
        return len(plot_order)

all_pngs = sorted(glob.glob("*.png"), key=plot_sort_key)
plot_b64 = {}
for png in all_pngs:
    with open(png, "rb") as f:
        plot_b64[png] = base64.b64encode(f.read()).decode()

plot_html = ""
if plot_b64:
    imgs = []
    for png, b64 in plot_b64.items():
        label = os.path.splitext(os.path.basename(png))[0].replace("_", " ").title()
        imgs.append(
            f'<figure><img src="data:image/png;base64,{b64}" '
            f'alt="{html.escape(label)}">'
            f'<figcaption>{html.escape(label)}</figcaption></figure>'
        )
    plot_html = "<h2>QC Plots</h2><div class='plots'>" + "".join(imgs) + "</div>"

now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
scope_warning = ""
if "${qc_scope}" == "provisional":
    scope_warning = "<p class='warn'>Sample-level QC is provisional because not all autosomes were analysed.</p>"

html_doc = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SNP Array QC Report — ${meta.id}</title>
<style>
body {{ font-family: Arial, sans-serif; max-width: 1400px; margin: auto; padding: 20px; color: #333; }}
h1 {{ color: #2c6fad; border-bottom: 2px solid #2c6fad; }}
h2 {{ color: #2c6fad; margin-top: 32px; }}
table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
th {{ background: #2c6fad; color: white; padding: 8px 12px; text-align: left; }}
td {{ border: 1px solid #ccc; padding: 7px 12px; }}
tr:nth-child(even) {{ background: #f9f9f9; }}
.info {{ background: #e3f2fd; padding: 12px; border-radius: 4px; margin: 10px 0; }}
.warn {{ color: #a33; }}
.plots {{ display: flex; flex-wrap: wrap; gap: 16px; }}
figure {{ margin: 0; flex: 0 0 calc(50% - 8px); box-sizing: border-box; }}
figure img {{ width: 100%; border: 1px solid #ddd; border-radius: 4px; }}
figcaption {{ text-align: center; font-size: 0.85em; color: #555; margin-top: 4px; }}
.scroll-table {{ max-height: 420px; overflow-y: auto; border: 1px solid #ccc; }}
.scroll-table table {{ margin-bottom: 0; }}
.scroll-table th {{ position: sticky; top: 0; z-index: 1; }}
td.yes {{ color: #b22; font-weight: bold; }}
td.no  {{ color: #555; }}
</style>
</head>
<body>
<h1>SNP Array QC Report</h1>
<p class="info"><strong>Dataset:</strong> ${meta.id} &nbsp;|&nbsp;
<strong>Chromosomes:</strong> ${params.chroms} &nbsp;|&nbsp;
<strong>Sample QC scope:</strong> ${qc_scope} &nbsp;|&nbsp;
<strong>Generated:</strong> {now}</p>
{scope_warning}

<h2>Summary</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Final samples</td><td>{final_samples}</td></tr>
<tr><td>Final variants</td><td>{final_variants}</td></tr>
<tr><td>Excluded samples</td><td>{excluded_sample_count}</td></tr>
<tr><td>Excluded variants</td><td>{excluded_variant_count}</td></tr>
</table>

<h2>QC Attrition</h2>
<table>
<tr><th>Phase</th><th>Step</th><th>Metric</th><th>Value</th></tr>
{''.join(summary_rows) if summary_rows else '<tr><td colspan="4">No summary files were produced.</td></tr>'}
</table>

<h2>Thresholds and Run Settings</h2>
<table>
<tr><th>Parameter</th><th>Value</th></tr>
{''.join(threshold_rows)}
</table>

{batch_html}

{per_sample_html}

{plot_html}
</body>
</html>
'''

with open("final_report.html", "w") as out:
    out.write(html_doc)

print(f"HTML report written with {len(plot_b64)} embedded plots")
PYEOF

    # ── PDF report ────────────────────────────────────────────────────────────
    # Assemble all QC plots into a multi-page PDF.
    # Uses base R (grid) + the png package (standard in most R installations).
    # Falls back gracefully if png package is absent.
    if command -v Rscript &>/dev/null; then
        Rscript - << 'RSCRIPT'
pngs <- sort(list.files(".", pattern="\\.png$", full.names=TRUE))

plot_order <- c("miss_het_scatter","heterozygosity_plot","sex_check_F_stat",
                "vmiss_plot","cc_miss_plot","maf_plot","hwe_plot","imputation_r2_plot",
                "relatedness_pi_hat","pca_scree","pca_plot")
key <- function(p) {
    n <- tools::file_path_sans_ext(basename(p))
    idx <- match(n, plot_order)
    ifelse(is.na(idx), length(plot_order) + 1L, idx)
}
pngs <- pngs[order(sapply(pngs, key))]

if (length(pngs) == 0) {
    message("No PNG plots found — skipping PDF generation")
    quit(save="no")
}

has_png_pkg <- requireNamespace("png", quietly=TRUE)
if (!has_png_pkg) {
    message("R package 'png' not available — skipping PDF generation")
    quit(save="no")
}

pdf("qc_report.pdf", width=10, height=8)

# Title page
grid::grid.newpage()
grid::grid.text("SNP Array QC Report",    x=0.5, y=0.65, gp=grid::gpar(fontsize=24, fontface="bold", col="#2c6fad"))
grid::grid.text("Dataset: ${meta.id}",    x=0.5, y=0.50, gp=grid::gpar(fontsize=14))
grid::grid.text(paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
                x=0.5, y=0.40, gp=grid::gpar(fontsize=11, col="grey50"))

# One plot per page
for (f in pngs) {
    img <- tryCatch(png::readPNG(f), error=function(e) NULL)
    if (is.null(img)) next
    grid::grid.newpage()
    grid::grid.raster(img, width=0.95, height=0.90, y=0.48)
    label <- gsub("_", " ", tools::file_path_sans_ext(basename(f)))
    grid::grid.text(label, x=0.5, y=0.97, gp=grid::gpar(fontsize=10, col="grey40"))
}

# Per-batch summary page
batch <- tryCatch(
    read.table("qc_per_batch.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE),
    error=function(e) NULL
)
if (!is.null(batch) && nrow(batch) > 0) {
    grid::grid.newpage()
    grid::grid.text("Per-batch QC Summary",
                    x=0.5, y=0.95,
                    gp=grid::gpar(fontsize=16, fontface="bold", col="#2c6fad"))
    hdrs <- c("Batch", "Samples", "Removed", "Pass", "Pass %")
    xs   <- c(0.05, 0.30, 0.45, 0.58, 0.72)
    y0   <- 0.87
    for (j in seq_along(hdrs))
        grid::grid.text(hdrs[j], x=xs[j], y=y0, just="left",
                        gp=grid::gpar(fontsize=11, fontface="bold"))
    grid::grid.lines(x=c(0.04, 0.88), y=rep(y0 - 0.025, 2),
                     gp=grid::gpar(col="grey60"))
    for (i in seq_len(nrow(batch))) {
        y <- y0 - 0.05 - (i - 1) * 0.055
        vals <- c(batch$BATCH[i], batch$N_SAMPLES[i], batch$N_REMOVED[i],
                  batch$N_PASS[i], paste0(batch$PASS_RATE_PCT[i], "%"))
        for (j in seq_along(vals))
            grid::grid.text(vals[j], x=xs[j], y=y, just="left",
                            gp=grid::gpar(fontsize=10))
    }
}

dev.off()
message(paste("PDF report written:", length(pngs), "plots"))
RSCRIPT
    fi
    """
}
