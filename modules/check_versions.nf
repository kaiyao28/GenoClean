/*
================================================================================
  MODULE: CHECK_VERSIONS (shared between SNP array and WGS/WES pipelines)
================================================================================
  Purpose:
    Print the version of every tool used by the pipeline at run start.
    Output is captured in a versions.txt file and displayed in the log.
    This makes runs reproducible and helps diagnose environment problems.

  Include from either pipeline:
    include { CHECK_VERSIONS } from '../modules/check_versions'

  The process checks each tool and silently skips any that are not installed
  rather than failing, because some tools are optional (R, VerifyBamID2).
================================================================================
*/

process CHECK_VERSIONS {
    label 'process_low'
    publishDir "${params.outdir}/logs", mode: params.publish_dir_mode

    output:
    path "versions.txt", emit: versions

    script:
    """
    # ── Collect tool versions ─────────────────────────────────────────────────
    {
        echo "=== Genetic QC Pipeline — Tool Versions ==="
        echo "Run date: \$(date)"
        echo ""

        # Nextflow version is available as an environment variable
        echo "nextflow: \${NXF_VER:-unknown}"

        # PLINK 1.9
        if command -v plink &>/dev/null; then
            echo "plink: \$(plink --version 2>&1 | head -1)"
        else
            echo "plink: NOT FOUND"
        fi

        # PLINK2
        if command -v plink2 &>/dev/null; then
            echo "plink2: \$(plink2 --version 2>&1 | head -1)"
        else
            echo "plink2: NOT FOUND"
        fi

        # samtools
        if command -v samtools &>/dev/null; then
            echo "samtools: \$(samtools --version 2>&1 | head -1)"
        else
            echo "samtools: NOT FOUND"
        fi

        # bcftools
        if command -v bcftools &>/dev/null; then
            echo "bcftools: \$(bcftools --version 2>&1 | head -1)"
        else
            echo "bcftools: NOT FOUND"
        fi

        # Picard
        if command -v picard &>/dev/null; then
            echo "picard: \$(picard --version 2>&1 | head -1)"
        else
            echo "picard: NOT FOUND"
        fi

        # FastQC
        if command -v fastqc &>/dev/null; then
            echo "fastqc: \$(fastqc --version 2>&1 | head -1)"
        else
            echo "fastqc: NOT FOUND"
        fi

        # mosdepth
        if command -v mosdepth &>/dev/null; then
            echo "mosdepth: \$(mosdepth --version 2>&1 | head -1)"
        else
            echo "mosdepth: NOT FOUND"
        fi

        # GATK
        if command -v gatk &>/dev/null; then
            echo "gatk: \$(gatk --version 2>&1 | head -1)"
        else
            echo "gatk: NOT FOUND"
        fi

        # VerifyBamID2 (optional)
        if command -v VerifyBamID &>/dev/null; then
            echo "verifybamid2: \$(VerifyBamID 2>&1 | grep -i version | head -1 || echo installed)"
        else
            echo "verifybamid2: NOT FOUND (will use GATK CalculateContamination fallback)"
        fi

        # tabix
        if command -v tabix &>/dev/null; then
            echo "tabix: \$(tabix --version 2>&1 | head -1)"
        else
            echo "tabix: NOT FOUND"
        fi

        # Python
        if command -v python3 &>/dev/null; then
            echo "python3: \$(python3 --version 2>&1)"
            python3 -c "import pandas, numpy, scipy; print(f'  pandas={pandas.__version__} numpy={numpy.__version__} scipy={scipy.__version__}')" 2>/dev/null || true
        else
            echo "python3: NOT FOUND"
        fi

        # R (optional)
        if command -v Rscript &>/dev/null; then
            echo "R: \$(Rscript --version 2>&1 | head -1)"
            Rscript -e "cat('  ggplot2=', as.character(packageVersion('ggplot2')), '\n')" 2>/dev/null || true
        else
            echo "R: NOT FOUND (QC plots will be skipped)"
        fi

        echo ""
        echo "=== End of version check ==="
    } > versions.txt

    # Also print to log so it's visible during the run
    cat versions.txt
    """
}
