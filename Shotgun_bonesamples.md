# Bone samples shotgun processing. Prepare barcodes for each samples in a barcodes.txt file
---
## 1. Demultiplex
`Barcodes txt file example`
```
AGCTTCATCT      GCATATGTGT      EBL003
GGATCTATCT      TGACTCGTGT      EBL024
AGCTTCATCT      ATGCGAGTGT      EBL028
TAGACGATCT      GCATATGTGT      M18918
GGATCTATCT      ATGCGAGTGT      M24018
AGCTTCATCT      TGACTCGTGT      M5653
```
`Slurm script for stacks demultiplexing`
```
#!/bin/bash -e
#SBATCH --cpus-per-task  16
#SBATCH --job-name       demult
#SBATCH --mem            50G
#SBATCH --time           24:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
module purge
module load Stacks/2.67-GCC-12.3.0
process_shortreads -1 forward_L2.fastq.gz -2 reverse_L2.fastq.gz \
                -i gzfastq \
                -b SRW_bone_barcodes.txt --index_index \
                -c -q

```
## 2. Remove barcodes and adapters and perform FASTQC analysis with TrimGalore
```
#!/bin/bash -e

#SBATCH --cpus-per-task  16
#SBATCH --job-name       qctrimg
#SBATCH --mem            20G
#SBATCH --time           24:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --array          1-6

module purge
module load TrimGalore/0.6.10-gimkl-2022a-Python-3.11.3-Perl-5.34.1

# Define the array of input files
declare -a FQ_FILES1=(
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/EBL003.1.fq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/EBL024.1.fq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/EBL028.1.fq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/M18918.1.fq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/M24018.1.fq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/M5653.1.fq.gz"
)

declare -a FQ_FILES2=(
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/EBL003.2.fq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/EBL024.2.fq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/EBL028.2.fq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/M18918.2.fq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/M24018.2.fq.gz"
"/nesi/nobackup/uoa02626/SRW_Sebastian/Shotgun/M5653.2.fq.gz"
)

declare -a OUTPUT_PREFIXES=(
  "EBL003"
  "EBL024"
  "EBL028"
  "M18918"
  "M24018"
  "M5653"
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
## 3. Mapping to ref genome or mitogenome to extract mitochondrial sequences or nuclear sequences.
`Slurm script for bwa`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  16
#SBATCH --job-name       demult
#SBATCH --mem            75G
#SBATCH --time           32:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --array          1-6

module purge
module load BWA/0.7.18-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0

REF_GENOME="/nesi/nobackup/uoa02626/SRW_Sebastian/Euau_mitogenome.fasta"

# Define the array of SFS files and their corresponding output prefixes
declare -a FQ_FILES1=(
  "trimmed_EBL003.1.fq.gz"
  "trimmed_EBL028.1.fq.gz"
  "trimmed_M24018.1.fq.gz"
  "trimmed_EBL024.1.fq.gz" 
  "trimmed_M18918.1.fq.gz"
  "trimmed_M5653.1.fq.gz"
)

declare -a FQ_FILES2=(
  "trimmed_EBL003.2.fq.gz"
  "trimmed_EBL028.2.fq.gz"
  "trimmed_M24018.2.fq.gz"
  "trimmed_EBL024.2.fq.gz"
  "trimmed_M18918.2.fq.gz"
  "trimmed_M5653.2.fq.gz"
)

declare -a OUTPUT_PREFIXES=(
  "mitoaln_Euau_EBL003"
  "mitoaln_Euau_EBL028"
  "mitoaln_Euau_M24018"
  "mitoaln_Euau_EBL024"
  "mitoaln_Euau_M18918"
  "mitoaln_Euau_M5653"
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
## 4. Generate consensus mitogenome sequence for each sample using SAMtools.
`Slurm script for SAMtools`
```
#!/bin/bash -e

#SBATCH --cpus-per-task  4
#SBATCH --job-name       mitocons
#SBATCH --mem            10
#SBATCH --time           4:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --array          1-6

module purge
module load SAMtools/1.19-GCC-12.3.0

declare -a BAM_FILES=(
  "mito_EBL003.bam"
  "mito_EBL028.bam"
  "mito_M24018.bam"
  "mito_EBL024.bam" 
  "mito_M18918.bam"
  "mito_M5653.bam"
)

declare -a OUTPUT_PREFIXES=(
  "mito_EBL003.fasta"
  "mito_EBL028.fasta"
  "mito_M24018.fasta"
  "mito_EBL024.fasta"
  "mito_M18918.fasta"
  "mito_M5653.fasta"
)

# Get the index for the current array job
INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Get the corresponding input files and output prefix
BAM=${BAM_FILES1[$INDEX]}
OUTPUT_PREFIX=${OUTPUT_PREFIXES[$INDEX]}

samtools consensus -f fasta -o $OUTPUT_PREFIX $REF $BAM
```
