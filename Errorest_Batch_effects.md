# Subset bam files to just autosomes.
---
## 1. FASTQC analysis, adapter and barcode removal using TrimGalore
`Slurm script for Trimgalore`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  1
#SBATCH --job-name       subset21
#SBATCH --mem            20G
#SBATCH --time           24:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --array          1-34

module purge
module load SAMtools/1.22-GCC-12.3.0

samples=(
BREA14.sorted.bam
Chile_SRW1.sorted.bam
Eau0113.sorted.bam
Eau0208.sorted.bam
Eau04EA02.sorted.bam
Eau08AI032.sorted.bam
Eau08AI165.sorted.bam
Eau09AI050.sorted.bam
Eau09AI086.sorted.bam
Eau_09_SA_F12.sorted.bam
Eau12SAF_04.sorted.bam
Eau16AF_22.sorted.bam
Eau_16S_AF_25.sorted.bam
Eau23WA09.sorted.bam
Eau23WA10.sorted.bam
Eau23WA23.sorted.bam
Eau23WA26.sorted.bam
Eau24_WA_001.sorted.bam
Eau24_WA_004.sorted.bam
Eau24_WA_007.sorted.bam
Eau24_WA_011.sorted.bam
Eau24_WA_025.sorted.bam
Eau94WA01.sorted.bam
Eau94WA04.sorted.bam
Eau94WA06.sorted.bam
Eau95WA14.sorted.bam
Euau_EBL003.sorted.bam
Euau_EBL024.sorted.bam
Euau_EBL028.sorted.bam
Euau_M18918.sorted.bam
Euau_M24018.sorted.bam
Euau_M5653.sorted.bam
Humpback.sorted.bam
MQ28.sorted.bam
)

bam=${samples[$SLURM_ARRAY_TASK_ID-1]}
prefix=${bam%.bam}

# Subsetting + clean header
samtools view -@ ${SLURM_CPUS_PER_TASK} -h "$bam" \
  HiC_scaffold_{1..21} | \
awk '
  /^@HD/ || /^@PG/ || /^@RG/ {print; next}
  /^@SQ/ {
    if ($2 ~ /^SN:HiC_scaffold_([1-9]$|1[0-9]$|20$|21$)/) print
    next
  }
  $3 ~ /^HiC_scaffold_([1-9]$|1[0-9]$|20$|21$)/ || $3=="*" {print}
' | \
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o "${prefix}.scaffolds1_21.bam"

# Index
samtools index "${prefix}.scaffolds1_21.bam"


```
## 2. Get consensus sequence of perfect sample and outgroup Humpback whale (humpback whale wgs reads aligned to SRW reference genome).
`Slurm script for ANGSD`
```
#!/bin/bash -e
#SBATCH --cpus-per-task  1
#SBATCH --job-name       est_error
#SBATCH --mem            30G
#SBATCH --time           1-00:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%A_%a.out
#SBATCH --error          %x_%A_%a.err
#SBATCH --hint           nomultithread

module purge
module load angsd/0.935-GCC-9.2.0

angsd \
  -i Eau23WA09.sorted.bam \
  -doCounts 1 \
  -doFasta 1 \
  -remove_bads 0 \
  -only_proper_pairs 0 \
  -uniqueOnly 0 \
  -minMapQ 30 \
  -minQ 30 \
  -nThreads ${SLURM_CPUS_PER_TASK} \
  -out Eau23WA09.consensus

```
## 3. Individual error estimation.
`Slurm script for ANGSD`
```
#!/bin/bash -e
#SBATCH --cpus-per-task  8
#SBATCH --job-name       est_error
#SBATCH --mem            100G
#SBATCH --time           1-00:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%A_%a.out
#SBATCH --error          %x_%A_%a.err
#SBATCH --hint           nomultithread

module purge
module load angsd/0.935-GCC-9.2.0

anc=/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/allBams/Humpback.consensus.fa.gz
bamlist=/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/allBams/filt21c_all_samples_bams.txt
ref=/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/allBams/Eau23WA09.consensus.fa.gz
scaffolds=/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/allBams/scaffolds.list

angsd -doAncError 2 \
   -anc ${anc} \
   -ref ${ref} \
   -nThreads 8 \
   -out wARG_alls_Eau23WA09_error.estimates \
   -bam $bamlist \
   -uniqueOnly 1 \
   -remove_bads 1 \
   -minMapQ 30 \
   -minQ 30 \
   -rf ${scaffolds}

```
