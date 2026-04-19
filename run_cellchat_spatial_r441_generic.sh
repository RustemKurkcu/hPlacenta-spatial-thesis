#!/bin/bash
#SBATCH --job-name=cellchat_r441
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=400G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=shan.kurkcu@uconn.edu
#SBATCH --output=/home/FCAM/skurkcu/projects/hPlacenta-spatial-thesis/output/logs/%x_%j.log
#SBATCH --error=/home/FCAM/skurkcu/projects/hPlacenta-spatial-thesis/output/logs/%x_%j.log

set -euo pipefail

PROJECT_ROOT="/home/FCAM/skurkcu/projects/hPlacenta-spatial-thesis"
SCRIPT_PATH="${SCRIPT_PATH:-${PROJECT_ROOT}/scripts/01_active_pipeline/03b_spatial_cellchat_allcells_min3_rdsfirst.R}"
CELLCHAT_OUT_ROOT="${CELLCHAT_OUT_ROOT:-output/cellchat/allcells_min3}"
CELLCHAT_INPUT_DIR="${CELLCHAT_INPUT_DIR:-output/objects}"
CELLCHAT_INPUT_BASENAME="${CELLCHAT_INPUT_BASENAME:-02_scored_misi_ido1}"

cd "${PROJECT_ROOT}"
mkdir -p output/logs output/objects output/reports "${CELLCHAT_OUT_ROOT}"

module purge
source /home/FCAM/skurkcu/miniconda3/etc/profile.d/conda.sh
conda activate /home/FCAM/skurkcu/miniconda3/envs/spatial_r441

export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export PATH="$CONDA_PREFIX/bin:$PATH"
export CELLCHAT_OUT_ROOT CELLCHAT_INPUT_DIR CELLCHAT_INPUT_BASENAME

echo "[$(date -Is)] PROJECT_ROOT=${PROJECT_ROOT}"
echo "[$(date -Is)] SCRIPT_PATH=${SCRIPT_PATH}"
echo "[$(date -Is)] CELLCHAT_OUT_ROOT=${CELLCHAT_OUT_ROOT}"
echo "[$(date -Is)] which R: $(which R)"
echo "[$(date -Is)] which Rscript: $(which Rscript)"
Rscript -e 'cat("R version:", R.version.string, "\n")'

Rscript -e '
pkgs <- c("Seurat","SeuratObject","dplyr","jsonlite","future","RANN","CellChat","SpatialCellChat")
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
print(ok)
if (!all(ok)) {
  cat("Missing packages:", paste(names(ok)[!ok], collapse=", "), "\n")
  quit(status = 1)
}
'

if [[ ! -f "${SCRIPT_PATH}" ]]; then
  echo "ERROR: Missing SCRIPT_PATH=${SCRIPT_PATH}" >&2
  exit 1
fi

echo "[$(date -Is)] Running ${SCRIPT_PATH}"
Rscript "${SCRIPT_PATH}"
echo "[$(date -Is)] Completed ${SCRIPT_PATH}"
