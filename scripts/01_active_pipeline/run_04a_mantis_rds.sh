#!/bin/bash
#SBATCH --job-name=04a_atlas_compute
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH --output=output/logs/slurm_04a_%j.log
#SBATCH --error=output/logs/slurm_04a_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p output/logs output/objects output/reports output/tables output/figures

module purge
module load R/4.2.2 cmake/3.20.2 gcc/10.2.0
export R_LIBS_USER="$HOME/rlibs"

# Force .rds object preference on cluster runs
export ATLAS_PREFER_RDS=1

# Optional knobs (override via sbatch --export=ALL,ATLAS_TARGET_CELLTYPE=...,ATLAS_MISI_LOW_Q=...,ATLAS_MISI_HIGH_Q=...)
: "${ATLAS_TARGET_CELLTYPE:=EVT}"
: "${ATLAS_MISI_LOW_Q:=0.33}"
: "${ATLAS_MISI_HIGH_Q:=0.67}"

echo "[$(date -Is)] Running 04a compute with RDS preference"
Rscript scripts/01_active_pipeline/04a_compute_spatial_cellchat_atlas.R
echo "[$(date -Is)] Completed 04a compute"
