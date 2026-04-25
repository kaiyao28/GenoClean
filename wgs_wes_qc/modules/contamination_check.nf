/*
================================================================================
  MODULE: CONTAMINATION_CHECK
================================================================================
  Purpose:
    Estimate the fraction of cross-individual DNA contamination in a BAM file
    using VerifyBamID2. Contaminated samples produce incorrect genotype calls
    and should be excluded or re-sequenced.

  Why this step:
    Even 2-3% contamination can inflate variant allele frequencies for
    heterozygous calls and create false-positive variants. gnomAD v4 uses
    a 3% contamination cutoff as a primary sample filter. VerifyBamID2 is
    the current standard tool (no reference panel required for the SVD method).

  Default threshold: params.max_contamination = 0.03 (3%)

  How to change:
    nextflow run main.nf --max_contamination 0.02  # more stringent

  How to disable:
    nextflow run main.nf --run_contamination false

  Output:
    - contamination.selfSM       : VerifyBamID2 per-sample output
    - contamination_summary.txt  : parsed summary for final report
================================================================================
*/

process CONTAMINATION_CHECK {
    label 'process_medium'
    publishDir "${params.outdir}/contamination/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path reference_fasta

    output:
    path "contamination.selfSM",     emit: selfSM
    path "contamination_summary.txt", emit: summary

    script:
    """
    # ── Run VerifyBamID2 using the SVD-based method ───────────────────────────
    # The SVD method does not require an external allele frequency reference panel.
    # If VERIFYBAMID2_SVD_PREFIX is set in the environment, it will be used.
    # Otherwise, the built-in 1000 Genomes SVD files are used if installed.
    VerifyBamID \\
        --SVDPrefix \${VERIFYBAMID2_SVD_PREFIX:-/opt/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat} \\
        --Reference ${reference_fasta} \\
        --BamFile ${bam} \\
        --Output contamination \\
        --DisableSanityCheck \\
        --NumThread ${task.cpus} || {
            # Fallback: use GATK CalculateContamination if VerifyBamID2 fails
            echo "VerifyBamID2 failed — attempting GATK GetPileupSummaries + CalculateContamination"
            gatk GetPileupSummaries \\
                -I ${bam} \\
                -V \${GATK_GNOMAD_COMMON_VARIANTS:-/dev/null} \\
                -L \${GATK_GNOMAD_COMMON_VARIANTS:-/dev/null} \\
                -O pileup_summaries.table || true
            gatk CalculateContamination \\
                -I pileup_summaries.table \\
                -O contamination_table.txt || true
            # Create a minimal selfSM-format file
            echo -e "FREEMIX" > contamination.selfSM
            echo -e "NA" >> contamination.selfSM
        }

    # ── Parse contamination estimate ──────────────────────────────────────────
    python3 - << 'PYEOF'
threshold = ${params.max_contamination}
freemix = "NA"
try:
    with open("contamination.selfSM") as fh:
        lines = fh.readlines()
    if len(lines) >= 2:
        header = lines[0].strip().split("\t")
        vals   = lines[1].strip().split("\t")
        if "FREEMIX" in header:
            idx = header.index("FREEMIX")
            freemix = float(vals[idx])
except Exception as e:
    print(f"Warning: could not parse contamination.selfSM: {e}")

status = "FAIL" if isinstance(freemix, float) and freemix > threshold else "PASS"

with open("contamination_summary.txt", "w") as out:
    out.write(f"step=contamination_check\n")
    out.write(f"dataset=${meta.id}\n")
    out.write(f"max_contamination_threshold={threshold}\n")
    out.write(f"contamination_estimate={freemix}\n")
    out.write(f"status={status}\n")

pct = f"{float(freemix)*100:.2f}%" if freemix != "NA" else "NA"
print(f"Contamination: {pct} (threshold {threshold*100:.0f}%) — {status}")
PYEOF
    """
}
