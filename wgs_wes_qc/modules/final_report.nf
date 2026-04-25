/*
================================================================================
  MODULE: FINAL_REPORT (WGS/WES)
================================================================================
  Purpose:
    Aggregate QC module summaries into one HTML report and TSV tables.
================================================================================
*/

process FINAL_REPORT {
    label 'process_report'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path summary_files

    output:
    path "wgs_wes_final_report.html", emit: report
    path "wgs_wes_qc_summary.tsv", emit: summary_table
    path "wgs_wes_thresholds.tsv", emit: thresholds

    script:
    """
    python3 - << 'PYEOF'
import datetime
import glob
import html

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

with open("wgs_wes_qc_summary.tsv", "w") as out:
    out.write("phase\\tstep\\tdataset\\tmetric\\tvalue\\n")
    for data in summaries:
        step = data.get("step", data["file"].replace(".txt", ""))
        if step in ("variant_calling_qc", "variant_filtering", "merge_chromosomes", "index_chrom_vcf", "index_input_vcf"):
            phase = "variant_level_qc"
        elif step == "input_check":
            phase = "input_check"
        else:
            phase = "sample_level_qc"
        dataset = data.get("dataset", "")
        for key, value in sorted(data.items()):
            if key in ("file", "step", "dataset"):
                continue
            out.write(f"{phase}\\t{step}\\t{dataset}\\t{key}\\t{value}\\n")

thresholds = [
    ("mode", "${params.mode}"),
    ("input_type", "${params.input_type}"),
    ("chroms", "${params.chroms}"),
    ("sample_qc_scope", "${params.sample_qc_scope}"),
    ("run_variant_qc", "${params.run_variant_qc}"),
    ("run_sample_qc", "${params.run_sample_qc}"),
    ("min_mean_depth_wgs", "${params.min_mean_depth_wgs}"),
    ("min_mean_depth_wes", "${params.min_mean_depth_wes}"),
    ("min_target_20x_fraction", "${params.min_target_20x_fraction}"),
    ("max_contamination", "${params.max_contamination}"),
    ("max_duplication_rate", "${params.max_duplication_rate}"),
    ("min_call_rate", "${params.min_call_rate}"),
    ("min_gq", "${params.min_gq}"),
    ("min_dp", "${params.min_dp}"),
    ("variant_qual", "${params.variant_qual}"),
    ("variant_filter_method", "${params.variant_filter_method}"),
    ("relatedness_pi_hat", "${params.relatedness_pi_hat}"),
    ("pca_outlier_sd", "${params.pca_outlier_sd}"),
]

with open("wgs_wes_thresholds.tsv", "w") as out:
    out.write("parameter\\tvalue\\n")
    for key, value in thresholds:
        out.write(f"{key}\\t{value}\\n")

rows = []
try:
    with open("wgs_wes_qc_summary.tsv") as fh:
        next(fh)
        for line in fh:
            phase, step, dataset, metric, value = line.rstrip("\\n").split("\\t", 4)
            rows.append(
                f"<tr><td>{html.escape(phase)}</td><td>{html.escape(step)}</td>"
                f"<td>{html.escape(dataset)}</td><td>{html.escape(metric)}</td>"
                f"<td>{html.escape(value)}</td></tr>"
            )
except Exception:
    pass

threshold_rows = [
    f"<tr><td>{html.escape(k)}</td><td>{html.escape(v)}</td></tr>"
    for k, v in thresholds
]

now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
html_doc = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>WGS/WES QC Report</title>
<style>
body {{ font-family: Arial, sans-serif; max-width: 1200px; margin: auto; padding: 20px; color: #333; }}
h1 {{ color: #2c6fad; border-bottom: 2px solid #2c6fad; }}
h2 {{ color: #2c6fad; margin-top: 32px; }}
table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
th {{ background: #2c6fad; color: white; padding: 8px 12px; text-align: left; }}
td {{ border: 1px solid #ccc; padding: 7px 12px; }}
tr:nth-child(even) {{ background: #f9f9f9; }}
.info {{ background: #e3f2fd; padding: 12px; border-radius: 4px; margin: 10px 0; }}
.warn {{ color: #a33; }}
</style>
</head>
<body>
<h1>WGS / WES QC Report</h1>
<p class="info"><strong>Mode:</strong> ${params.mode.toUpperCase()} &nbsp;|&nbsp;
<strong>Input type:</strong> ${params.input_type} &nbsp;|&nbsp;
<strong>Chromosomes:</strong> ${params.chroms} &nbsp;|&nbsp;
<strong>Generated:</strong> {now}</p>
<p class="warn">Sample-level QC should be treated as final only when all autosomes were analysed.</p>

<h2>QC Summary</h2>
<table>
<tr><th>Phase</th><th>Step</th><th>Dataset</th><th>Metric</th><th>Value</th></tr>
{''.join(rows) if rows else '<tr><td colspan="5">No summary files were produced.</td></tr>'}
</table>

<h2>Thresholds And Run Settings</h2>
<table>
<tr><th>Parameter</th><th>Value</th></tr>
{''.join(threshold_rows)}
</table>
</body>
</html>
'''

with open("wgs_wes_final_report.html", "w") as out:
    out.write(html_doc)
PYEOF
    """
}
