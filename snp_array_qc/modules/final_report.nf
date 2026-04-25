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
    val qc_scope

    output:
    path "final_report.html", emit: report
    path "qc_attrition_table.tsv", emit: attrition
    path "qc_thresholds.tsv", emit: thresholds

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
    ("hwe_p", "${params.hwe_p}"),
    ("maf", "${params.maf}"),
    ("heterozygosity_sd", "${params.heterozygosity_sd}"),
    ("relatedness_pi_hat", "${params.relatedness_pi_hat}"),
    ("sex_check_f_lower_female", "${params.sex_check_f_lower_female}"),
    ("sex_check_f_upper_male", "${params.sex_check_f_upper_male}"),
    ("pca_outlier_sd", "${params.pca_outlier_sd}"),
]

with open("qc_thresholds.tsv", "w") as out:
    out.write("parameter\\tvalue\\n")
    for key, value in thresholds:
        out.write(f"{key}\\t{value}\\n")

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

now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
scope_warning = ""
if "${qc_scope}" == "provisional":
    scope_warning = "<p class='warn'>Sample-level QC is provisional because not all autosomes were analysed.</p>"

html_doc = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SNP Array QC Report - ${meta.id}</title>
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

<h2>QC Steps</h2>
<table>
<tr><th>Phase</th><th>Step</th><th>Metric</th><th>Value</th></tr>
{''.join(summary_rows) if summary_rows else '<tr><td colspan="4">No summary files were produced.</td></tr>'}
</table>

<h2>Thresholds And Run Settings</h2>
<table>
<tr><th>Parameter</th><th>Value</th></tr>
{''.join(threshold_rows)}
</table>
</body>
</html>
'''

with open("final_report.html", "w") as out:
    out.write(html_doc)
PYEOF
    """
}
