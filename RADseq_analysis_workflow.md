# Scripts for population structure, kinship, and genetic diversity analysis using ddRADseq data.
---
## 1. Quality trimming using trimmomatic
`Slurm script for NGSadmix`
```
#!/bin/bash -e
#SBATCH --cpus-per-task  8
#SBATCH --job-name       qctrimg
#SBATCH --mem            8G
#SBATCH --time           04:00:00
#SBATCH --account        uoa02626
#SBATCH --output         logs/%x_%A_%a.out
#SBATCH --error          logs/%x_%A_%a.err
#SBATCH --array          1-346

# trim_ddrad_array.sl
# Trims all demultiplexed ddRAD individual fastq pairs with Trim Galore.
# Resources reduced relative to the WGS version (8 cpus/8G/4h vs 16 cpus/20G/24h)
# since per-individual ddRAD read files are far smaller than WGS lane data.
#
# Requires 3 list files with matching line order (line N in each file = same sample):
#   r1s_allfqs_list.txt   - full path to each sample's R1 fastq.gz
#   r2s_allfqs_list.txt   - full path to each sample's R2 fastq.gz
#   base_ids.txt          - output prefix (sample ID) for each sample
#
# IMPORTANT: verify all three lists are sorted/ordered identically before running -
# see the alignment check below, which aborts the whole job if lengths mismatch.

module purge
module load TrimGalore/0.6.10-gimkl-2022a-Python-3.11.3-Perl-5.34.1

LIST_DIR=/nesi/nobackup/uoa02626/FKW_Sebastian
R1_LIST="${LIST_DIR}/r1s_allfqs_list.txt"
R2_LIST="${LIST_DIR}/r2s_allfqs_list.txt"
ID_LIST="${LIST_DIR}/base_ids.txt"

OUT_DIR=/nesi/nobackup/uoa02626/FKW_Sebastian/trimmed_fqs
mkdir -p "${OUT_DIR}" logs

# Load the three lists into arrays
mapfile -t FQ_FILES1 < "${R1_LIST}"
mapfile -t FQ_FILES2 < "${R2_LIST}"
mapfile -t OUTPUT_PREFIXES < "${ID_LIST}"

# Sanity check: all three lists must be the same length or indices won't line up
if [[ ${#FQ_FILES1[@]} -ne ${#FQ_FILES2[@]} || ${#FQ_FILES1[@]} -ne ${#OUTPUT_PREFIXES[@]} ]]; then
    echo "ERROR: list length mismatch - R1=${#FQ_FILES1[@]} R2=${#FQ_FILES2[@]} IDs=${#OUTPUT_PREFIXES[@]}" >&2
    exit 1
fi

INDEX=$((SLURM_ARRAY_TASK_ID - 1))

FQ1=${FQ_FILES1[$INDEX]}
FQ2=${FQ_FILES2[$INDEX]}
OUTPUT_PREFIX=${OUTPUT_PREFIXES[$INDEX]}

echo "=== Trimming ${OUTPUT_PREFIX} ==="
echo "R1: ${FQ1}"
echo "R2: ${FQ2}"

trim_galore -q 20 -j 8 --fastqc --paired \
    -o "${OUT_DIR}" --basename "${OUTPUT_PREFIX}" \
    "${FQ1}" "${FQ2}"

echo "=== Done ${OUTPUT_PREFIX} ==="
```
## 4. Admixture analysis using NGSadmix
`Slurm script for NGSadmix`
```
#!/bin/bash -e
#SBATCH --cpus-per-task  1
#SBATCH --job-name       admixngs
#SBATCH --mem            50G
#SBATCH --time           1-00:00:00
#SBATCH --account        uoa02626
#SBATCH --output         %x_%j.out
#SBATCH --error          %x_%j.err
#SBATCH --hint           nomultithread

module purge
module load angsd/0.935-GCC-9.2.0

# Define variables
BEAGLE_FILE="/nesi/nobackup/uoa02626/SRW_Sebastian/RADseq/2_RAD_SRW_highc_angsd.beagle.gz"
OUT_PREFIX="/nesi/nobackup/uoa02626/SRW_Sebastian/RADseq/RADseq_ngsadmix"
THREADS=${SLURM_CPUS_PER_TASK}
MIN_MAF=0.05                  # Minimum MAF
MIS_TOL=0.8                  # Tolerance for high-quality genotypes
SEEDS=(35000)    # Seeds for reproducibility
MAX_ITER=2000                # Maximum number of EM iterations

# Loop through K values
for K in {2..8}; do
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
