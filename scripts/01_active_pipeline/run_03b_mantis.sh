#!/bin/bash
#SBATCH --job-name=03b_cellchat_full
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=400G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH --output=output/logs/slurm_03b_%j.log
#SBATCH --error=output/logs/slurm_03b_%j.log

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
mkdir -p output/logs output/objects output/reports

module purge
module load R/4.2.2 cmake/3.20.2 gcc/10.2.0

# CRITICAL: Point R to the personal package library
export R_LIBS_USER="$HOME/rlibs"

TARGET_SCRIPT="scripts/01_active_pipeline/03b_spatial_cellchat_all_pathways.R"
if [[ ! -f "${TARGET_SCRIPT}" ]]; then
  TARGET_SCRIPT="scripts/01_active_pipeline/03b_spatial_cellchat_allcells.R"
fi

if [[ ! -f "${TARGET_SCRIPT}" ]]; then
  echo "ERROR: Could not find 03b CellChat script. Checked all_pathways and allcells variants." >&2
  exit 1
fi

echo "[$(date -Is)] Running ${TARGET_SCRIPT}"
Rscript "${TARGET_SCRIPT}"
echo "[$(date -Is)] Completed ${TARGET_SCRIPT}"
