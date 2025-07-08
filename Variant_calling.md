# Raw sequences preprocessing.
---
## 1. FASTQC analysis, adapter and barcode removal using TrimGalore
`Slurm script for Trimgalore`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  16
#SBATCH --job-name       qctrimg
#SBATCH --mem            20G
#SBATCH --time           24:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --array          1-9

module purge
module load TrimGalore/0.6.10-gimkl-2022a-Python-3.11.3-Perl-5.34.1

REF_GENOME="/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/RWref_HiC.fasta"

# Define the array of input files
declare -a FQ_FILES1=(
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_1-Eau09AI050/sample_1_combined_R1.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_2-Eau08AI032/sample_2_combined_R1.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_3-Eau08AI165/sample_3_combined_R1.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_4-Eau09AI086/sample_4_combined_R1.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_5-Eau08AI006/sample_5_combined_R1.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_6-Eau94WA01/sample_6_combined_R1.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_7-Eau94WA04/sample_7_combined_R1.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_8-Eau94WA06/sample_8_combined_R1.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_9-Eau95WA14/sample_9_combined_R1.fastq.gz"
)

declare -a FQ_FILES2=(
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_1-Eau09AI050/sample_1_combined_R2.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_2-Eau08AI032/sample_2_combined_R2.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_3-Eau08AI165/sample_3_combined_R2.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_4-Eau09AI086/sample_4_combined_R2.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_5-Eau08AI006/sample_5_combined_R2.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_6-Eau94WA01/sample_6_combined_R2.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_7-Eau94WA04/sample_7_combined_R2.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_8-Eau94WA06/sample_8_combined_R2.fastq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/SRWwholegen/allSRW_WGS/Sample_9-Eau95WA14/sample_9_combined_R2.fastq.gz"
)

declare -a OUTPUT_PREFIXES=(
  "Eau09AI050"
  "Eau08AI032"
  "Eau08AI165"
  "Eau09AI086"
  "Eau08AI006"
  "Eau94WA01"
  "Eau94WA04"
  "Eau94WA06"
  "Eau95WA14"
)

# Get the index for the current array job
INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Get the corresponding input files and output prefix
FQ1=${FQ_FILES1[$INDEX]}
FQ2=${FQ_FILES2[$INDEX]}
OUTPUT_PREFIX=${OUTPUT_PREFIXES[$INDEX]}

# Ensure output directory exists

# Run Trim Galore
trim_galore -q 20 -j 16 --fastqc --paired \
    -o ./ --basename ${OUTPUT_PREFIX} \
    ${FQ1} ${FQ2}

```
## 2. Mapping to SRW reference genome using BWA.
`Slurm script for BWA`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  16
#SBATCH --job-name       alnbwa
#SBATCH --mem            75G
#SBATCH --time           32:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --array          1-9

module purge
module load BWA/0.7.18-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0

REF_GENOME="/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/RWref_HiC.fasta"

# Define the array of SFS files and their corresponding output prefixes
declare -a FQ_FILES1=(
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau08AI006_val_1.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau08AI032_val_1.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau08AI165_val_1.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau09AI050_val_1.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau09AI086_val_1.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau94WA01_val_1.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau94WA04_val_1.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau94WA06_val_1.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau95WA14_val_1.fq.gz"
)

declare -a FQ_FILES2=(
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau08AI006_val_2.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau08AI032_val_2.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau08AI165_val_2.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau09AI050_val_2.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau09AI086_val_2.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau94WA01_val_2.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau94WA04_val_2.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau94WA06_val_2.fq.gz"
  "/nesi/nobackup/uoa02626/SRW_Sebastian/WGS_SRW/Trimming/Eau95WA14_val_2.fq.gz"
)

declare -a OUTPUT_PREFIXES=(
  "Eau08AI032"
  "Eau08AI165"
  "Eau09AI050"
  "Eau09AI086"
  "Eau94WA01"
  "Eau94WA04"
  "Eau94WA06"
  "Eau95WA14"
)

# Get the index for the current array job
INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Get the corresponding input files and output prefix
FQ1=${FQ_FILES1[$INDEX]}
FQ2=${FQ_FILES2[$INDEX]}
OUTPUT_PREFIX=${OUTPUT_PREFIXES[$INDEX]}

# Run BWA mem and process the output with SAMtools
bwa mem -M -t 16 -R "@RG\tID:${OUTPUT_PREFIX}\tSM:${OUTPUT_PREFIX}" "$REF_GENOME" "$FQ1" "$FQ2" | \
samtools view -F 4 -Sb - | \
samtools sort -o "${OUTPUT_PREFIX}.sorted.bam" - 

# Index the sorted BAM file
samtools index -c "${OUTPUT_PREFIX}.sorted.bam"
```
## 3. Remove duplicates using GATK MarkDuplicates.
`Slurm script for GATK`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  4
#SBATCH --job-name       dedup
#SBATCH --mem            30G
#SBATCH --time           24:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --array          1-14

module purge
module load GATK/4.6.1.0-gimkl-2022a
module load SAMtools/1.19-GCC-12.3.0

# Define the array of BAM files and their corresponding output prefixes
declare -a IN_BAM=(
"Eau08AI032.sorted.bam"
"Eau09AI086.sorted.bam"
"Eau94WA06.sorted.bam"
"Euau_EBL024.sorted.bam"
"Euau_M24018.sorted.bam"
"Eau08AI165.sorted.bam"
"Eau94WA01.sorted.bam"
"Eau95WA14.sorted.bam"
"Euau_EBL028.sorted.bam"
"Euau_M5653.sorted.bam"
"Eau09AI050.sorted.bam"
"Eau94WA04.sorted.bam"
"Euau_EBL003.sorted.bam"
"Euau_M18918.sorted.bam"
)

declare -a OUTPUT_PREFIXES=(
"Eau08AI032"
"Eau09AI086"
"Eau94WA06"
"Euau_EBL024"
"Euau_M24018"
"Eau08AI165"
"Eau94WA01"
"Eau95WA14"
"Euau_EBL028"
"Euau_M5653"
"Eau09AI050"
"Eau94WA04"
"Euau_EBL003"
"Euau_M18918"
)

# Get the index for the current array job
INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Get the corresponding input files and output prefix
IB=${IN_BAM[$INDEX]}
OUTPUT_PREFIX=${OUTPUT_PREFIXES[$INDEX]}

# Run BWA mem and process the output with SAMtools
gatk --java-options "-Xmx20g" MarkDuplicates \
    -I "${IB}" \
    -O "${OUTPUT_PREFIX}_dedup.bam" \
    -M "${OUTPUT_PREFIX}_metrics.txt" \
    --MAX_RECORDS_IN_RAM 5000 \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --VALIDATION_STRINGENCY SILENT

samtools index -c "${OUTPUT_PREFIX}_dedup.bam"
```
## 4. Admixture analysis with NGSadmix from genotype likelihoods dataset. Tested K values from 2 to 10.
`Slurm script for NGSadmix`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  16
#SBATCH --job-name       admixngs
#SBATCH --mem            50G
#SBATCH --time           24:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load angsd/0.935-GCC-9.2.0

# Define variables
BEAGLE_FILE="nrlt_alldol_angsd_gl.beagle.gz"
OUT_PREFIX="nrlt_alldol_ngsAdmix_run"
THREADS=${SLURM_CPUS_PER_TASK}
MIN_MAF=0.1                  # Minimum MAF
MIS_TOL=0.8                  # Tolerance for high-quality genotypes
SEEDS=(11045 45690 59321)    # Seeds for reproducibility
MAX_ITER=2000                # Maximum number of EM iterations

# Loop through K values
for K in {2..10}; do
    echo "Running ngsAdmix for K=${K}"

    # Loop through seeds
    for SEED in "${SEEDS[@]}"; do
        echo "Using seed=${SEED}"

        # Run ngsAdmix
./NGSadmix -likes $BEAGLE_FILE \
                 -K $K \
                 -outfiles "${OUT_PREFIX}_K${K}_seed${SEED}" \
                 -seed $SEED \
                 -minMaf $MIN_MAF \
                 -misTol $MIS_TOL \
                 -P $THREADS \
                 -maxiter $MAX_ITER
    done
done
```
## 5. PCAs using genotype likelihoods datasets. We performed a PCA using all Māui and Hector's individuals and a PCA using just Hector's dataset.
`Slurm script for PCAngsd`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  16
#SBATCH --job-name       angsdpca
#SBATCH --mem            25G
#SBATCH --time           4:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load angsd/0.935-GCC-9.2.0

conda activate myenv #Python environment where pcangsd was installed

pcangsd -b nrlt_alldol_angsd_gl.beagle.gz \
   --maf 0.1 \
   --tree \
   --inbreed_samples \
   -o nrlt_pca \
   -t ${SLURM_CPUS_PER_TASK}
```
## 6. Admixture result evaluation using EvalAdmix.
`Slurm script for evalAdmix`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  16
#SBATCH --job-name       admixeval
#SBATCH --mem            20G
#SBATCH --time           8:00:00
#SBATCH --account        uoo02423
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load angsd/0.935-GCC-9.2.0

./evalAdmix -beagle nrlt_alldol_angsd_gl.beagle.gz \
     -fname nrlt_alldol_ngsAdmix_run_K6.fopt.gz \
     -qname nrlt_alldol_ngsAdmix_run_K6.qopt \
     -P ${SLURM_CPUS_PER_TASK} \
     -o nrlt_ngsAdmix_alldol.corres.txt
```
