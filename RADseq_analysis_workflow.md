# Scripts for population structure, kinship, and genetic diversity analysis using ddRADseq data.
---
## 1. Admixture analysis using NGSadmix
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
