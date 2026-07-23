# Runs of homozygosity and Genome wide heteroygosity analyzis from a VCF file after filtering and preprocessing with VCFtools and BCFtools.
---
## 1. VCF file rehead, annotate and filter
```
module load BCFtools/1.22-GCC-12.3.0
module load VCFtools/0.1.17-GCC-12.3.0-Perl-5.38.2

#Reformat vcf file
bcftools view \
    -e 'GT="."' \
    -m2 -M2 \
    -v snps \
    -Oz \
    -o ALL_srw_outfile.fil.vcf.gz \
    NEW_fil_SRW_highc.vcf.gz

bcftools annotate \
    --header-line '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">' \
    -Oz \
    -o ALL_fixed_header.vcf.gz \
    ALL_srw_outfile.fil.vcf.gz

#Filter MAF, missing data and keep only biallelic snps
vcftools --vcf ALL_fixed_header.vcf.gz \
   --recode-INFO-all \
   --maf 0.01 \
   --max-alleles 2 \
   --min-alleles 2 \
   --minDP 5 \
   --max-missing 0.9 \
   --out ALLW_filtered \
   --recode \
   --remove-indels
```
## 2. ROH analysis using bcftools roh (best when variable coverage across samples)
`Slurm script for bcftools roh`
```                                                                                             
#!/bin/bash -e
#SBATCH --cpus-per-task  1
#SBATCH --job-name       rohbcftools
#SBATCH --mem            30G
#SBATCH --time           24:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load BCFtools/1.22-GCC-12.3.0

# Generate list of samples, exclude related individuals for calculating allele frequencies
# in bcftools roh command later

FIXED="ALL_fixed_header.vcf.gz"
ROH_OUT="ALL_srw_roh"

# 1. Compute allele frequencies
echo "[1] Computing allele frequencies..."
bcftools +fill-tags "${FIXED}" -Oz -o ALL_fixed_header.AF.vcf.gz -- -t AF
bcftools index --tbi ALL_fixed_header.AF.vcf.gz

# 2. Run bcftools roh
echo "[2] Running bcftools roh..."
bcftools roh \
    --AF-tag AF \
    --rec-rate 1e-8 \
    -G30 \
    -Or \
    -o "${ROH_OUT}_full.txt" \
    ALL_fixed_header.AF.vcf.gz

# 3. Extract ROH segments only (RG records)
echo "[3] Extracting ROH segments..."
grep "^RG" "${ROH_OUT}_full.txt" > "${ROH_OUT}_segments.txt"
echo "Done. ROH segments: ${ROH_OUT}_segments.txt"
echo "Total segments: $(wc -l < ${ROH_OUT}_segments.txt)"
```
## 3. ROH analysis using PLINK
`Slurm script for PLINK`
```
#!/bin/bash -e
#SBATCH --cpus-per-task  4
#SBATCH --job-name       plinkroh
#SBATCH --mem            60G
#SBATCH --time           24:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load PLINK/1.09b6.16

plink --vcf ALL_fixed_header.vcf.gz \
      --homozyg \
      --homozyg-snp 30 \
      --homozyg-window-snp 30 \
      --homozyg-window-het 1 \
      --homozyg-window-missing 5 \
      --homozyg-density 50 \
      --homozyg-gap 500 \
      --homozyg-kb 100 \
      --allow-extra-chr \
      --double-id \
      --out ALL_srw_plink_ROH
```
## 4. Run GWH with script "gw_slidwin.py", for this script we need a file with chromosome sizes of the reference genome. Run an array with the number of chromosomes, it runs each chromosome independently.
Get chromosome sizes from reference genome .fai file
`cut -f1,2 RWref_HiC.fasta.fai > chrom_sizes.txt`
`Slurm script to run custom python script`
```
#!/bin/bash -e
#SBATCH --cpus-per-task  1
#SBATCH --job-name       gwh_py
#SBATCH --mem            10G
#SBATCH --time           3:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread
#SBATCH --array          1-21

module load Python/3.11.6-foss-2023a

# Define the path to the chromosome_lengths.txt file
CHROM_LENGTHS="chrom_sizes.txt"

# Define other constants
VCF_FILE="ALL_fixed_header.vcf.gz"
WINDOW_SIZE=2000000
STEP_SIZE=200000

# Read the chromosome number for the current array task, change this to the start of the chromosome id on the ref genome, "Chromosome_"
chrom=HiC_scaffold_$SLURM_ARRAY_TASK_ID

# Run the Python script for the chromosome
python ./gw_slidwin.py \
    "$VCF_FILE" \
    "$CHROM_LENGTHS" \
    "$WINDOW_SIZE" \
    "$STEP_SIZE" \
    "$chrom"
```

