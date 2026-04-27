# Setup Guide

This guide explains how to set up the Genetic QC Nextflow pipeline on a laptop, workstation, or HPC cluster.

The key idea is:

```text
Nextflow runs on your computer.
Docker or Apptainer provides the bioinformatics tools.
```

So you need both:

```text
1. Nextflow on the host machine
2. A container engine:
   - Docker for laptops/workstations
   - Apptainer or Singularity for HPC clusters
```

## What Each Tool Does

| Tool | Why it is needed |
|------|------------------|
| Git | Downloads the pipeline repository. |
| Java | Required by Nextflow. Use Java 17 or newer. |
| Nextflow | Runs the workflow and launches each QC step. |
| Docker | Runs the prebuilt container image on laptops/workstations. |
| Apptainer/Singularity | Runs the container image on HPC clusters where Docker is usually not allowed. |

The Docker image contains tools such as PLINK, PLINK2, bcftools, samtools, GATK, FastQC, mosdepth, Python, and R.

## Quick Decision Guide

| Where are you running? | Recommended setup |
|------------------------|-------------------|
| Windows laptop | Docker Desktop on Windows + Nextflow inside WSL Ubuntu |
| Linux workstation | Nextflow + Docker |
| macOS workstation | Nextflow + Docker Desktop |
| HPC cluster | Nextflow + Apptainer/Singularity + scheduler profile |

## Windows Setup

Recommended Windows setup:

```text
PowerShell:
  install/start Docker Desktop

WSL Ubuntu:
  install Java
  install Nextflow
  clone and run the pipeline
```

PowerShell is not Bash. Commands such as `curl -s ... | bash`, `export PATH=...`, and `bash test_env.sh` are Linux/Bash commands. Use them inside WSL Ubuntu, not plain PowerShell.

### 1. Install Docker Desktop

Install Docker Desktop:

```text
https://www.docker.com/products/docker-desktop/
```

Open Docker Desktop and wait until it says the engine is running.

In PowerShell, check:

```powershell
docker version
docker info
```

Pull the published image:

```powershell
docker pull ghcr.io/kaiyao28/genetic-qc:1.0
```

If this fails with `dockerDesktopLinuxEngine`, Docker Desktop is not running or the Linux/WSL backend is not enabled.

### 2. Install WSL Ubuntu

In PowerShell:

```powershell
wsl.exe --install Ubuntu
```

Restart Windows if prompted. Then open Ubuntu from the Start Menu.

### 3. Install Java, Git, and Nextflow inside Ubuntu

Inside Ubuntu:

```bash
sudo apt update
sudo apt install -y openjdk-17-jre curl git
```

Install Nextflow:

```bash
curl -s https://get.nextflow.io | bash
mkdir -p ~/bin
mv nextflow ~/bin/
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Check:

```bash
java -version
nextflow -version
```

### 4. Enable Docker access from WSL

In Docker Desktop:

```text
Settings -> Resources -> WSL Integration -> enable Ubuntu
```

Back in Ubuntu:

```bash
docker version
docker pull ghcr.io/kaiyao28/genetic-qc:1.0
```

### 5. Clone and run the test workflow

Inside Ubuntu:

```bash
git clone https://github.com/kaiyao28/GeneticQC.git
cd GeneticQC
```

Run the smoke tests:

```bash
bash test_data/run_smoke_tests.sh                    # both pipelines
bash test_data/run_smoke_tests.sh --test snp_array   # SNP array only
bash test_data/run_smoke_tests.sh --test wgs_wes     # WGS/WES only
```

Each test runs the selected workflow on synthetic toy data in `test_data/` and writes an HTML report to `results/`.

## Linux or macOS Setup

Install Java 17 or newer, Git, Docker, and Nextflow.

Check:

```bash
java -version
git --version
docker version
nextflow -version
```

Install Nextflow if needed:

```bash
curl -s https://get.nextflow.io | bash
mkdir -p ~/bin
mv nextflow ~/bin/
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

Clone the pipeline:

```bash
git clone https://github.com/kaiyao28/GeneticQC.git
cd GeneticQC
```

Pull the Docker image:

```bash
docker pull ghcr.io/kaiyao28/genetic-qc:1.0
```

Run the smoke tests:

```bash
bash test_data/run_smoke_tests.sh                    # both pipelines
bash test_data/run_smoke_tests.sh --test snp_array   # SNP array only
bash test_data/run_smoke_tests.sh --test wgs_wes     # WGS/WES only
```

## HPC Cluster Setup

Most clusters do not allow Docker. Use Apptainer or Singularity instead.

Typical requirements:

```text
Java
Nextflow
Apptainer or Singularity
Scheduler profile, for example SLURM
```

Check available modules:

```bash
module avail java
module avail nextflow
module avail apptainer
module avail singularity
```

Load modules. Exact names vary by cluster:

```bash
module load java
module load nextflow
module load apptainer
```

If Nextflow is not provided as a module, install it in your home directory:

```bash
curl -s https://get.nextflow.io | bash
mkdir -p ~/bin
mv nextflow ~/bin/
export PATH="$HOME/bin:$PATH"
```

Create or pull a SIF image on shared storage:

```bash
mkdir -p containers
apptainer pull containers/genetic-qc.sif docker://ghcr.io/kaiyao28/genetic-qc:1.0
```

If your cluster uses `singularity` instead of `apptainer`:

```bash
mkdir -p containers
singularity pull containers/genetic-qc.sif docker://ghcr.io/kaiyao28/genetic-qc:1.0
```

Make sure `conf/singularity.config` points to the SIF file. The default is:

```groovy
process {
    container = "${projectDir}/containers/genetic-qc.sif"
}
```

If the SIF is elsewhere, use an absolute shared path:

```groovy
process {
    container = "/shared/containers/genetic-qc.sif"
}
```

Run with the cluster profile, for example SLURM plus Singularity:

```bash
nextflow run wgs_wes_qc/main.nf \
  --input_type vcf \
  --samplesheet samplesheet.csv \
  --reference_fasta /shared/reference/GRCh38.fa \
  --mode wgs \
  --chroms 1-22 \
  --outdir /shared/results/wgs_wes_qc \
  -profile slurm,singularity \
  -resume
```

Ask your cluster support team:

```text
Which scheduler is used: SLURM, LSF, PBS?
Which module loads Java, Nextflow, and Apptainer/Singularity?
Where should Nextflow work directories go?
Which storage path is shared across compute nodes?
Are internet pulls from compute/login nodes allowed?
```

## HPC Without Containers (Manual Tool Installation)

Use this path when the cluster has no Docker, no Apptainer/Singularity, and no
conda/mamba — but you can download files on the login node.

### Requirements

You need these already available (most HPC clusters have them):

```text
Java 17+    — for Nextflow, GATK, Picard, FastQC
gcc + make  — to compile samtools, bcftools, htslib from source
wget or curl
```

Check before starting:

```bash
java -version
gcc --version
module avail java    # if java is missing, load it first
module avail gcc     # if gcc is missing, load it first
```

### Step 1 — Download all tools

Run the setup script on the login node. It downloads pre-compiled binaries
where possible and compiles from source otherwise. No root access needed.

```bash
# Default: installs to ~/genetic_qc_tools
bash setup_hpc_manual.sh

# Custom location (recommended for shared clusters):
bash setup_hpc_manual.sh --dir /shared/software/genetic_qc_tools

# Skip samtools/bcftools compilation (if they are already available via module load):
bash setup_hpc_manual.sh --skip-compile
```

The script downloads: Nextflow, PLINK 1.9, PLINK 2, samtools, bcftools, htslib,
GATK 4, Picard, FastQC, mosdepth, and VerifyBamID2. All executables are placed
under `$TOOL_DIR/bin/`.

### Step 2 — Point the pipeline at those tools

Set an environment variable so Nextflow knows where the tools live:

```bash
export GENETIC_QC_TOOL_DIR=~/genetic_qc_tools   # or your custom path
```

Add this to `~/.bashrc` so it persists across sessions:

```bash
echo 'export GENETIC_QC_TOOL_DIR=~/genetic_qc_tools' >> ~/.bashrc
```

### Step 3 — Run the pipeline with the manual_paths profile

```bash
# Local execution
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --outdir results/snp_array_qc \
  -profile manual_paths

# SLURM cluster
nextflow run wgs_wes_qc/main.nf \
  --input_type vcf \
  --samplesheet samplesheet.csv \
  --reference_fasta /shared/reference/GRCh38.fa \
  --mode wgs \
  --outdir results/wgs_wes_qc \
  -profile slurm,manual_paths \
  -resume

# LSF cluster
nextflow run wgs_wes_qc/main.nf \
  --input_type bam \
  --samplesheet samplesheet.csv \
  --reference_fasta /shared/reference/GRCh38.fa \
  --mode wes \
  --target_intervals /shared/reference/exome_targets.bed \
  --outdir results/wgs_wes_qc \
  -profile lsf,manual_paths \
  -resume
```

Alternatively, pass `--tool_dir` directly without setting the environment variable:

```bash
nextflow run snp_array_qc/main.nf \
  --bfile data/raw/genotypes \
  --tool_dir /shared/software/genetic_qc_tools \
  -profile slurm,manual_paths
```

### Troubleshooting manual installs

**samtools/bcftools compilation fails**

Load gcc first, then re-run the script:

```bash
module load gcc
bash setup_hpc_manual.sh
```

Or ask the cluster sysadmin to load samtools and bcftools as modules, then use
`--skip-compile` and add the module to your job submission script:

```bash
module load samtools bcftools
bash setup_hpc_manual.sh --skip-compile
```

**VerifyBamID2 download fails**

This is non-fatal. The contamination module falls back to GATK
CalculateContamination automatically. You will see a warning in the log.

**A specific tool is already installed as a cluster module**

You can load it in your shell before running Nextflow and it will be available:

```bash
module load gatk samtools
export GENETIC_QC_TOOL_DIR=~/genetic_qc_tools
nextflow run wgs_wes_qc/main.nf -profile slurm,manual_paths ...
```

Tools on `PATH` when Nextflow starts are inherited by compute node processes.

## Testing The Environment

### Step 1 — Verify tools

Docker:

```bash
bash test_env.sh docker
```

For a different Docker image:

```bash
GENETIC_QC_DOCKER_IMAGE=my-image:tag bash test_env.sh docker
```

Singularity/Apptainer:

```bash
bash test_env.sh singularity
```

Manual install (`manual_paths` profile):

```bash
export PATH=$GENETIC_QC_TOOL_DIR/bin:$PATH
plink --version && plink2 --version && samtools --version | head -1
bcftools --version | head -1 && gatk --version && mosdepth --version
```

On Windows PowerShell without WSL/Git Bash, test the Docker image directly:

```powershell
docker run --rm ghcr.io/kaiyao28/genetic-qc:1.0 plink --version
docker run --rm ghcr.io/kaiyao28/genetic-qc:1.0 bcftools --version
docker run --rm ghcr.io/kaiyao28/genetic-qc:1.0 gatk --version
```

### Step 2 — Generate test data and run smoke tests

The synthetic test data is tracked in the repository. If you need to regenerate it
(e.g. after resetting the repo or changing the data):

```bash
python3 test_data/generate_test_data.py
rm -f test_data/snp_array/toy.bed test_data/snp_array/toy.bim test_data/snp_array/toy.fam
```

Run both pipelines:

```bash
bash test_data/run_smoke_tests.sh                          # Docker (default)
bash test_data/run_smoke_tests.sh --profile singularity    # Singularity
bash test_data/run_smoke_tests.sh --profile manual_paths   # no container
```

Run one pipeline at a time:

```bash
bash test_data/run_smoke_tests.sh --profile manual_paths --test snp_array
bash test_data/run_smoke_tests.sh --profile manual_paths --test wgs_wes
```

The test data is designed so QC filters have something to remove:
SNP array has monomorphic variants (fail MAF) and high-missingness variants (fail callrate).
WGS/WES VCF has variants with failing site-level flags and samples with low genotype quality.

## Common Problems

### `nextflow: command not found`

Nextflow is not installed or not on `PATH`.

Check:

```bash
nextflow -version
```

If missing, install Nextflow and add it to `PATH`.

### PowerShell `curl -s ... | bash` fails

PowerShell is not Bash. Use WSL Ubuntu, Git Bash, or Linux/macOS.

### Docker says `dockerDesktopLinuxEngine`

Docker Desktop is not running, or WSL2/Linux backend is not enabled.

Open Docker Desktop and check:

```powershell
docker version
docker info
```

### Docker fails with `input/output error` while pulling

This usually means Docker Desktop's internal Linux storage is full, partially
corrupted, or contains a broken layer from a previous failed pull.

First restart Docker Desktop. Then try:

```bash
docker image rm ghcr.io/kaiyao28/genetic-qc:1.0
docker builder prune
docker system prune
```

Then re-run:

```bash
bash test_data/run_smoke_tests.sh
# or run one pipeline at a time:
bash test_data/run_smoke_tests.sh --test snp_array
bash test_data/run_smoke_tests.sh --test wgs_wes
```

If it still fails, increase Docker Desktop disk space:

```text
Docker Desktop -> Settings -> Resources -> Advanced -> Disk image size
```

As a last resort:

```text
Docker Desktop -> Troubleshoot -> Clean / Purge data
```

This removes local Docker images and containers, but not your Git repository.

### `docker pull ghcr.io/kaiyao28/genetic-qc:1.0` says denied

The GHCR package is private or the image has not been published. Maintainers should check GitHub Actions and package visibility.

### Cluster does not allow Docker

Use Apptainer/Singularity and the `singularity` profile instead of Docker.
