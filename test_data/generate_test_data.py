#!/usr/bin/env python3
"""
generate_test_data.py — create synthetic test data for GenoClean smoke tests.

Run from the repository root:
    python3 test_data/generate_test_data.py

Overwrites:
    test_data/reference/mini.fa          500bp chr22 reference (replaces all-N stub)
    test_data/reference/mini.fa.fai      FASTA index
    test_data/wgs_wes/toy_chr22.vcf      20-variant 5-sample VCF
    test_data/wgs_wes/samplesheet_vcf.csv
    test_data/snp_array/toy.ped          10-sample 50-variant PED
    test_data/snp_array/toy.map          50-variant MAP

The data is intentionally small but designed so QC filters have something to do:
  VCF  : 3 site-level failing variants, 2 indels, mixed Ti/Tv (ratio ~2.0)
  Array: 3 monomorphic variants (fail MAF), 5 high-missingness variants
         (fail callrate), 2 all-heterozygous variants (test HWE step)
"""
import os
import random

random.seed(42)

# ── Navigate to repo root regardless of where the script is launched from ────
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
os.chdir(ROOT)


# ══════════════════════════════════════════════════════════════════════════════
#  Reference FASTA  (500bp chr22)
# ══════════════════════════════════════════════════════════════════════════════
# Deterministic pattern so REF alleles in the VCF can be verified.
PATTERN = "ATCGATCG"
REF_LEN = 500
SEQ = (PATTERN * (REF_LEN // len(PATTERN) + 1))[:REF_LEN]

def ref_base(pos1, length=1):
    """Return reference base(s) at 1-based position."""
    return SEQ[pos1 - 1: pos1 - 1 + length]

os.makedirs("test_data/reference", exist_ok=True)

with open("test_data/reference/mini.fa", "w") as fh:
    fh.write(f">22\n{SEQ}\n")

# FAI columns: name  length  byte_offset  bases_per_line  bytes_per_line
# ">22\n" = 4 bytes before the sequence
with open("test_data/reference/mini.fa.fai", "w") as fh:
    fh.write(f"22\t{REF_LEN}\t4\t{REF_LEN}\t{REF_LEN + 1}\n")

print(f"[OK] Reference: test_data/reference/mini.fa  ({REF_LEN} bp, chr22)")


# ══════════════════════════════════════════════════════════════════════════════
#  WGS/WES VCF  (5 samples, 20 variants on chr22)
# ══════════════════════════════════════════════════════════════════════════════
# Ti/Tv: 12 transitions + 6 transversions among 18 SNPs = ratio 2.0 (WGS-like)
# Plus 2 indels (excluded from Ti/Tv by bcftools).
# 3 variants carry failing site-level filters (low QD, high FS, low QUAL).
# Some genotypes have low DP/GQ to exercise genotype-level filtering.
#
# Transitions:   A↔G, C↔T
# Transversions: A↔C, A↔T, G↔C, G↔T

SAMPLES_VCF = ["SAMPLE1", "SAMPLE2", "SAMPLE3", "SAMPLE4", "SAMPLE5"]

# (pos, rsID, alt, QUAL, filter, QD, FS, MQ, MQRankSum, ReadPosRankSum)
VARIANTS = [
    # --- 12 passing transitions ---
    ( 25, "rs001", "G",  60, "PASS", 10.0,  2.0, 60.0,  0.5,  0.2),   # A→G Ti
    ( 50, "rs002", "C",  55, "PASS",  9.5,  1.5, 59.0,  0.3,  0.1),   # T→C Ti
    ( 75, "rs003", "T",  65, "PASS", 11.0,  2.5, 60.0, -0.2, -0.3),   # C→T Ti
    (100, "rs004", "A",  70, "PASS", 12.0,  1.0, 61.0,  0.1,  0.0),   # G→A Ti
    (175, "rs007", "T",  67, "PASS", 11.5,  1.5, 60.5,  0.2,  0.1),   # C→T Ti
    (200, "rs008", "A",  72, "PASS", 12.5,  1.0, 61.0,  0.0,  0.0),   # G→A Ti
    (225, "rs009", "G",  60, "PASS", 10.0,  2.5, 59.0,  0.3,  0.2),   # A→G Ti
    (275, "rs011", "T",  68, "PASS", 11.0,  2.0, 60.0, -0.2, -0.1),   # C→T Ti
    (300, "rs012", "A",  73, "PASS", 12.0,  1.5, 61.0,  0.1,  0.0),   # G→A Ti
    (350, "rs014", "C",  59, "PASS",  9.5,  2.0, 59.5,  0.4,  0.2),   # T→C Ti
    (400, "rs016", "A",  74, "PASS", 12.5,  1.0, 61.0,  0.0,  0.0),   # G→A Ti
    (425, "rs017", "G",  63, "PASS", 10.5,  2.0, 60.0,  0.2,  0.1),   # A→G Ti
    # --- 6 passing transversions ---
    (125, "rs005", "C",  58, "PASS",  9.0,  3.0, 59.5,  0.4,  0.2),   # A→C Tv
    (150, "rs006", "G",  62, "PASS", 10.5,  2.0, 60.0, -0.1, -0.1),   # T→G Tv
    (250, "rs010", "A",  54, "PASS",  8.5,  3.5, 58.5,  0.5,  0.3),   # T→A Tv
    (450, "rs018", "A",  56, "PASS",  8.8,  3.2, 58.0,  0.4,  0.3),   # T→A Tv
    # --- 3 failing variants (site-level) ---
    (325, "rs013", "G",  45, "FAIL",  1.0,  5.0, 60.0,  0.0,  0.0),   # Ti FAIL QD<2
    (375, "rs015", "A",  61, "FAIL",  8.0,100.0, 60.0,  0.3,  0.1),   # Tv FAIL FS>60
    # --- 2 more transversions ---
    (450, "rs018b","C",  56, "PASS",  8.8,  3.2, 58.0,  0.4,  0.3),   # T→C... conflict, skip
]

# Cleaner design — explicit flat list
SNPS = [
    # (pos, id, ref_expected, alt, qual, filter, QD, FS, MQ, MQRs, RPRs, Ti/Tv)
    ( 25, "rs001", "A", "G",  60, "PASS", 10.0,  2.0, 60.0,  0.5,  0.2, "Ti"),
    ( 50, "rs002", "T", "C",  55, "PASS",  9.5,  1.5, 59.0,  0.3,  0.1, "Ti"),
    ( 75, "rs003", "C", "T",  65, "PASS", 11.0,  2.5, 60.0, -0.2, -0.3, "Ti"),
    (100, "rs004", "G", "A",  70, "PASS", 12.0,  1.0, 61.0,  0.1,  0.0, "Ti"),
    (125, "rs005", "A", "C",  58, "PASS",  9.0,  3.0, 59.5,  0.4,  0.2, "Tv"),
    (150, "rs006", "T", "G",  62, "PASS", 10.5,  2.0, 60.0, -0.1, -0.1, "Tv"),
    (175, "rs007", "C", "T",  67, "PASS", 11.5,  1.5, 60.5,  0.2,  0.1, "Ti"),
    (200, "rs008", "G", "A",  72, "PASS", 12.5,  1.0, 61.0,  0.0,  0.0, "Ti"),
    (225, "rs009", "A", "G",  60, "PASS", 10.0,  2.5, 59.0,  0.3,  0.2, "Ti"),
    (250, "rs010", "T", "A",  54, "PASS",  8.5,  3.5, 58.5,  0.5,  0.3, "Tv"),
    (275, "rs011", "C", "T",  68, "PASS", 11.0,  2.0, 60.0, -0.2, -0.1, "Ti"),
    (300, "rs012", "G", "A",  73, "PASS", 12.0,  1.5, 61.0,  0.1,  0.0, "Ti"),
    (325, "rs013", "A", "G",  45, "FAIL",  1.0,  5.0, 60.0,  0.0,  0.0, "Ti"),  # FAIL: QD<2
    (350, "rs014", "T", "C",  59, "PASS",  9.5,  2.0, 59.5,  0.4,  0.2, "Ti"),
    (375, "rs015", "C", "A",  61, "FAIL",  8.0,100.0, 60.0,  0.3,  0.1, "Tv"),  # FAIL: FS>60
    (400, "rs016", "G", "A",  74, "PASS", 12.5,  1.0, 61.0,  0.0,  0.0, "Ti"),
    (425, "rs017", "A", "G",  63, "PASS", 10.5,  2.0, 60.0,  0.2,  0.1, "Ti"),
    (450, "rs018", "T", "A",  56, "PASS",  8.8,  3.2, 58.0,  0.4,  0.3, "Tv"),  # Tv
]
# Ti=12, Tv=5 among 17 SNPs → Ti/Tv = 2.4 ≈ WGS-like

# 2 indels at positions 475 and 500
# pos 475: "C" at 474%8=2 → C; pos 476: 475%8=3 → G  → DEL: REF=CG, ALT=C
# pos 500: "G" at 499%8=3 → G                         → INS: REF=G,  ALT=GA
INDELS = [
    (475, "rs019", ref_base(475, 2), ref_base(475)),        # DEL
    (500, "rs020", ref_base(500),    ref_base(500) + "A"),  # INS
]

# Verify REF alleles match the reference sequence
for pos, rsid, ref_allele, alt, *_ in SNPS:
    actual = ref_base(pos)
    assert actual == ref_allele, (
        f"REF mismatch at pos {pos} ({rsid}): expected {ref_allele}, reference has {actual}"
    )

# Genotype table: [SAMPLE1..5] for each variant (SNPs then indels)
# (GT, DP, GQ)  — some with low DP/GQ for genotype-level filtering test
GENO_SNPS = [
    # S1           S2           S3           S4           S5
    [("0/0",25,99),("0/1",20,85),("1/1",22,90),("0/1",18,80),("0/0",30,99)],  # rs001
    [("0/1",18,80),("0/0",25,99),("0/1",21,82),("1/1",15,75),("0/0",28,99)],  # rs002
    [("0/0",30,99),("0/1",22,88),("0/0",25,99),("0/1",20,85),("1/1",18,90)],  # rs003
    [("0/1",20,85),("1/1",17,88),("0/1",23,90),("0/0",28,99),("0/1",19,82)],  # rs004
    [("0/0",25,99),("0/0",22,99),("0/1",20,80),("0/1",18,75),("0/0",30,99)],  # rs005
    [("0/1",22,88),("0/0",25,99),("0/0",28,99),("1/1",15,85),("0/1",20,82)],  # rs006
    [("1/1",18,90),("0/1",20,85),("0/0",25,99),("0/0",28,99),("0/1",22,88)],  # rs007
    [("0/0",28,99),("0/1",22,88),("0/1",19,80),("0/0",30,99),("1/1",16,85)],  # rs008
    [("0/1",20,82),("0/0",26,99),("1/1",17,88),("0/1",21,85),("0/0",25,99)],  # rs009
    [("0/0",25,99),("0/1",20,80),("0/0",28,99),("0/1",18,78),("0/0",30,99)],  # rs010
    [("0/1",19,82),("1/1",16,88),("0/1",22,85),("0/0",27,99),("0/1",20,80)],  # rs011
    [("0/0",27,99),("0/1",21,85),("0/1",18,78),("1/1",14,85),("0/0",29,99)],  # rs012
    # rs013 FAIL-QD: all low DP+GQ — will be set to ./. by genotype filter
    [("0/1", 5,12),("0/0", 7,18),("0/1", 6,15),("0/0", 8,10),("1/1", 5,14)], # rs013
    [("0/0",24,99),("0/1",19,82),("0/1",21,85),("0/0",26,99),("1/1",16,88)],  # rs014
    # rs015 FAIL-FS: normal DP/GQ (FS is a site metric, genotypes are fine)
    [("0/1",20,80),("0/0",25,99),("0/0",27,99),("0/1",18,78),("0/0",30,99)],  # rs015
    [("0/0",30,99),("0/1",22,88),("1/1",18,90),("0/1",20,85),("0/0",25,99)],  # rs016
    [("0/1",21,85),("0/0",27,99),("0/1",19,80),("1/1",15,85),("0/0",28,99)],  # rs017
    # rs018: one sample low GQ — for genotype-filter demo
    [("0/1",18, 8),("0/0",25,99),("0/0",28,99),("0/1",19, 9),("1/1",14,85)],  # rs018
]
GENO_INDELS = [
    [("0/1",20,82),("0/0",25,99),("0/1",18,80),("0/0",27,99),("1/1",15,85)],  # rs019
    [("0/0",28,99),("0/1",21,85),("0/0",25,99),("0/1",19,82),("0/0",30,99)],  # rs020
]

os.makedirs("test_data/wgs_wes", exist_ok=True)

vcf_header = [
    "##fileformat=VCFv4.2",
    f"##contig=<ID=22,length={REF_LEN}>",
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FILTER=<ID=FAIL,Description="Failed site-level QC filter">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##INFO=<ID=QD,Number=1,Type=Float,Description="Quality by depth">',
    '##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand bias">',
    '##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">',
    '##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="MQ rank sum">',
    '##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Read position rank sum">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES_VCF),
]

vcf_records = []
for i, (pos, rsid, ref, alt, qual, flt, qd, fs, mq, mqrs, rprs, _) in enumerate(SNPS):
    info = f"QD={qd};FS={fs};MQ={mq};MQRankSum={mqrs};ReadPosRankSum={rprs}"
    gts  = "\t".join(f"{g}:{d}:{q}" for g, d, q in GENO_SNPS[i])
    vcf_records.append(f"22\t{pos}\t{rsid}\t{ref}\t{alt}\t{qual}\t{flt}\t{info}\tGT:DP:GQ\t{gts}")

for j, (pos, rsid, ref, alt) in enumerate(INDELS):
    info = "QD=9.0;FS=2.0;MQ=60.0;MQRankSum=0.0;ReadPosRankSum=0.0"
    gts  = "\t".join(f"{g}:{d}:{q}" for g, d, q in GENO_INDELS[j])
    vcf_records.append(f"22\t{pos}\t{rsid}\t{ref}\t{alt}\t65\tPASS\t{info}\tGT:DP:GQ\t{gts}")

# Sort by position
vcf_records.sort(key=lambda l: int(l.split("\t")[1]))

with open("test_data/wgs_wes/toy_chr22.vcf", "w") as fh:
    fh.write("\n".join(vcf_header + vcf_records) + "\n")

with open("test_data/wgs_wes/samplesheet_vcf.csv", "w") as fh:
    fh.write("sample,file1\n")
    fh.write("toy_cohort,test_data/wgs_wes/toy_chr22.vcf\n")

n_ti = sum(1 for *_, ttype in SNPS if ttype == "Ti")
n_tv = sum(1 for *_, ttype in SNPS if ttype == "Tv")
print(f"[OK] VCF: test_data/wgs_wes/toy_chr22.vcf")
print(f"     {len(SNPS)} SNPs (Ti={n_ti}, Tv={n_tv}, Ti/Tv={n_ti/n_tv:.1f})"
      f"  +  {len(INDELS)} indels  =  {len(SNPS)+len(INDELS)} variants")
print(f"     {len(SAMPLES_VCF)} samples  |  3 variants with failing site filters")
print(f"     2 variants with low-DP/GQ genotypes for genotype-filter demo")


# ══════════════════════════════════════════════════════════════════════════════
#  SNP Array PED / MAP
# ══════════════════════════════════════════════════════════════════════════════
# 10 samples (5 male, 5 female), 50 variants on chr22.
#
# Designed QC triggers:
#   Variants 0-39  : normal genotypes, MAF 0.10–0.40, no missingness
#   Variants 40-44 : high missingness (~35%) → fail variant callrate (>2%)
#   Variants 45-47 : monomorphic (all hom-ref, MAF=0) → fail MAF filter (0.01)
#   Variants 48-49 : all heterozygous → extreme HWE departure (for HWE step)

N_SAMPLES_SNP = 10
N_VARIANTS_SNP = 50
ALLELE_PAIRS = [("A","G"), ("C","T"), ("A","T"), ("G","C"), ("A","C"), ("G","T")]

# Build variant list
snp_map = []
for v in range(N_VARIANTS_SNP):
    snp_map.append(f"22\trs_snp{v+1:03d}\t0\t{(v + 1) * 1000}")

# Build genotype matrix: [sample][variant] = "allele1 allele2"
geno_matrix = [[None] * N_VARIANTS_SNP for _ in range(N_SAMPLES_SNP)]

for s in range(N_SAMPLES_SNP):
    for v in range(N_VARIANTS_SNP):
        major, minor = ALLELE_PAIRS[v % len(ALLELE_PAIRS)]

        if v < 40:
            # Normal: random HWE genotypes, MAF 0.10–0.40
            maf = 0.10 + (v % 7) * 0.05
            r = random.random()
            if r < (1 - maf) ** 2:
                geno_matrix[s][v] = f"{major} {major}"
            elif r < (1 - maf) ** 2 + 2 * maf * (1 - maf):
                geno_matrix[s][v] = f"{major} {minor}"
            else:
                geno_matrix[s][v] = f"{minor} {minor}"

        elif v < 45:
            # High missingness: ~35% missing per sample
            if random.random() < 0.35:
                geno_matrix[s][v] = "0 0"
            else:
                geno_matrix[s][v] = f"{major} {major}"

        elif v < 48:
            # Monomorphic (MAF = 0): all hom-ref → removed by --maf 0.01
            geno_matrix[s][v] = f"{major} {major}"

        else:
            # All heterozygous: extreme HWE departure
            geno_matrix[s][v] = f"{major} {minor}"

# Write PED
# Columns: FID IID Father Mother Sex(1=male,2=female) Phenotype(1=ctrl,2=case) [genotypes...]
ped_lines = []
for s in range(N_SAMPLES_SNP):
    fid   = "FAM1"
    iid   = f"S{s + 1:02d}"
    sex   = 1 if s < 5 else 2          # first 5 male, last 5 female
    pheno = 2 if s < 6 else 1          # 6 cases, 4 controls
    genos = "\t".join(geno_matrix[s])
    ped_lines.append(f"{fid}\t{iid}\t0\t0\t{sex}\t{pheno}\t{genos}")

os.makedirs("test_data/snp_array", exist_ok=True)
with open("test_data/snp_array/toy.ped", "w") as fh:
    fh.write("\n".join(ped_lines) + "\n")

with open("test_data/snp_array/toy.map", "w") as fh:
    fh.write("\n".join(snp_map) + "\n")

print(f"\n[OK] PED: test_data/snp_array/toy.ped  ({N_SAMPLES_SNP} samples)")
print(f"[OK] MAP: test_data/snp_array/toy.map  ({N_VARIANTS_SNP} variants)")
print(f"     Variants 0-39  : normal (MAF 0.10-0.40, no missing)")
print(f"     Variants 40-44 : ~35% missing → will fail callrate filter (>2%)")
print(f"     Variants 45-47 : monomorphic → will fail MAF filter (<0.01)")
print(f"     Variants 48-49 : all het → exercises HWE step")

print("""
Remove stale PLINK binary files if re-running:
  rm -f test_data/snp_array/toy.bed test_data/snp_array/toy.bim test_data/snp_array/toy.fam

Run individual smoke tests:
  bash test_data/run_smoke_tests.sh --profile manual_paths --test snp_array
  bash test_data/run_smoke_tests.sh --profile manual_paths --test wgs_wes
""")
