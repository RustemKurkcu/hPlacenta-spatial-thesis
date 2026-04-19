#!/bin/bash
#SBATCH --job-name=04b_atlas_plot
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH --output=output/logs/slurm_04b_%j.log
#SBATCH --error=output/logs/slurm_04b_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p output/logs output/objects output/reports output/tables output/figures

module purge
module load R/4.2.2 cmake/3.20.2 gcc/10.2.0
export R_LIBS_USER="$HOME/rlibs"

# Force .rds object preference on cluster runs
export ATLAS_PREFER_RDS=1

echo "[$(date -Is)] Running 04b plotting with RDS preference"
Rscript scripts/01_active_pipeline/04b_plot_spatial_cellchat_atlas.R
echo "[$(date -Is)] Completed 04b plotting"
