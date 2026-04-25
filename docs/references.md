# References

## SNP Array QC

**Anderson CA, Pettersson FH, Clarke GM, Cardon LR, Morris AP, Zondervan KT (2010)**
Data quality control in genetic case-control association studies.
*Nature Protocols*, 5(9): 1564–1573.
https://doi.org/10.1038/nprot.2010.116

The primary reference for GWAS QC. Describes sample and variant missingness filters,
sex checks, heterozygosity, relatedness (IBD), and HWE filtering. Thresholds used in
this pipeline are taken directly from Table 1 of this paper.

---

**Turner SD, Armstrong LL, Bradford Y, Carlson CS, Crawford DC, Crenshaw AT, et al. (2011)**
Quality control procedures for genome-wide association studies.
*Current Protocols in Human Genetics*, Chapter 1, Unit 1.19.
https://doi.org/10.1002/0471142905.hg0119s68

Comprehensive GWAS QC tutorial. Covers marker-level and sample-level QC, population
stratification (PCA), and imputation QC.

---

**Marees AT, de Kluiver H, Stringer S, Vorspan F, Curis E, Marie-Claire C, Derks EM (2018)**
A tutorial on conducting genome-wide association studies: Quality control and statistical
analysis.
*International Journal of Methods in Psychiatric Research*, 27(2): e1608.
https://doi.org/10.1002/mpr.1608

Step-by-step tutorial implementing the Anderson et al. protocol with PLINK commands.
Includes heterozygosity and PCA steps with example code.

---

**PLINK 1.9 documentation**
Purcell S, Chang C. PLINK 1.9. http://www.cog-genomics.org/plink/1.9/

Reference for all PLINK commands used in the SNP array pipeline:
- `--mind`: sample missingness
- `--geno`: variant missingness
- `--hwe`: Hardy-Weinberg equilibrium filter
- `--maf`: minor allele frequency filter
- `--het`: heterozygosity
- `--genome`: pairwise IBD
- `--check-sex`: sex check
- `--pca`: principal component analysis
- `--indep-pairwise`: LD pruning

---

## WGS / WES QC

**McKenna A, Hanna M, Banks E, et al. (2010)**
The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA
sequencing data.
*Genome Research*, 20(9): 1297–1303.
https://doi.org/10.1101/gr.107524.110

Original GATK paper. The variant annotation thresholds used in this pipeline (QD, FS,
MQ, MQRankSum, ReadPosRankSum) are documented in the GATK framework.

---

**GATK Best Practices for germline short variant discovery**
https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices

The authoritative pipeline reference for WGS/WES variant calling. Covers BQSR, variant
calling with HaplotypeCaller, joint genotyping, VQSR, and hard filtering.

---

**GATK Hard-filtering germline short variants**
https://gatk.broadinstitute.org/hc/en-us/articles/360035890471

Specific guidance on hard-filter thresholds for SNPs and indels when VQSR is not feasible.
The thresholds implemented as defaults in this pipeline come directly from this document:
- SNPs: QD < 2.0, FS > 60, MQ < 40, MQRankSum < −12.5, ReadPosRankSum < −8
- Indels: QD < 2.0, FS > 200, ReadPosRankSum < −20

---

**gnomAD v4.0 QC methods**
Karczewski KJ, et al. (2023)
gnomAD v4.0 release notes.
https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/

Describes sample-level and variant-level QC applied to the gnomAD v4 release, including:
- Contamination cutoff (FREEMIX > 3%)
- Coverage-based sex inference
- Population PCA and outlier detection
- Per-sample variant count QC

---

**VerifyBamID2**
Zhang F, Flickinger M, Taliun SAG, et al. (2020)
Ancestry-agnostic estimation of DNA sample contamination from sequence reads.
*Genome Research*, 30(2): 185–194.
https://doi.org/10.1101/gr.246934.118

Tool used for cross-sample contamination estimation. The SVD-based method does not
require a population-specific allele frequency reference panel.

---

**Picard tools**
Broad Institute. Picard toolkit. http://broadinstitute.github.io/picard/

Used for:
- `MarkDuplicates`: PCR duplicate detection and marking
- `CollectAlignmentSummaryMetrics`: mapping rate and chimeric read metrics
- `CollectInsertSizeMetrics`: insert size distribution

---

**mosdepth**
Pedersen BS, Quinlan AR (2018)
mosdepth: quick coverage calculation for genomes and exomes.
*Bioinformatics*, 34(5): 867–868.
https://doi.org/10.1093/bioinformatics/btx699

Fast, memory-efficient coverage calculation tool used for WGS/WES coverage QC.

---

**bcftools**
Danecek P, Bonfield JK, Liddle J, et al. (2021)
Twelve years of SAMtools and BCFtools.
*GigaScience*, 10(2): giab008.
https://doi.org/10.1093/gigascience/giab008

Used throughout the WGS/WES pipeline for VCF statistics, filtering, and format conversion.

---

## Population structure / PCA

**Price AL, Patterson NJ, Plenge RM, et al. (2006)**
Principal components analysis corrects for stratification in genome-wide association studies.
*Nature Genetics*, 38(8): 904–909.
https://doi.org/10.1038/ng1847

Foundational paper establishing PCA as the standard method for controlling population
stratification in GWAS.

---

**1000 Genomes Project Consortium (2015)**
A global reference for human genetic variation.
*Nature*, 526(7571): 68–74.
https://doi.org/10.1038/nature15393

The 1000 Genomes Phase 3 dataset is commonly used as a reference panel for ancestry PCA
projection. The pipeline supports merging study data with a user-supplied reference panel.
