#!/bin/bash
set -euo pipefail

PROJECT_ROOT="/home/FCAM/skurkcu/projects/hPlacenta-spatial-thesis"
RUNNER="${PROJECT_ROOT}/run_cellchat_spatial_r441_generic.sh"

STD_SCRIPT="${PROJECT_ROOT}/scripts/01_active_pipeline/03b_spatial_cellchat_allcells_min3_rdsfirst.R"
ENH_SCRIPT="${PROJECT_ROOT}/scripts/01_active_pipeline/03b_spatial_cellchat_allcells_enhanced_toxicswitch.R"

std_job=$(sbatch --job-name=cellchat_std \
  --export=ALL,SCRIPT_PATH="${STD_SCRIPT}",CELLCHAT_OUT_ROOT="output/cellchat/allcells_min3",CELLCHAT_INPUT_DIR="output/objects",CELLCHAT_INPUT_BASENAME="02_scored_misi_ido1" \
  "${RUNNER}" | awk '{print $4}')

enh_job=$(sbatch --job-name=cellchat_enh \
  --export=ALL,SCRIPT_PATH="${ENH_SCRIPT}",CELLCHAT_OUT_ROOT="output/cellchat/enhanced_toxicswitch_min3",CELLCHAT_INPUT_DIR="output/objects",CELLCHAT_INPUT_BASENAME="02_scored_misi_ido1" \
  "${RUNNER}" | awk '{print $4}')

echo "Submitted standard job: ${std_job}"
echo "Submitted enhanced job: ${enh_job}"
echo "Watch queue with: squeue -u $USER"
