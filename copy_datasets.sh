#!/bin/bash
# copy_datasets.sh
# Copy all datasets used by bootstrap_ghf_component.m into a clean, timestamped folder

# base destination folder
BASE="datasets_copy"
STAMP=$(date +"%Y-%m-%d_%H-%M-%S")
DEST="${BASE}_${STAMP}"

# create destination
mkdir -p "$DEST"

# list of dataset files (from your MATLAB script)
FILES=(
  "data_ipics/GHF_Hazzard_interp.mat"
  "data_ipics/GHF_Martos_interp.mat"
  "data_ipics/GHF_Shen_interp.mat"
  "data_ipics/GHF_Stal_interp.mat"
  "data_ipics/GHF_Losing_interp.mat"
  "data_ipics/GHF_An_interp.mat"
  "data_ipics/GHF_FoxMaule_interp.mat"
  "data_ipics/UNC_Hazzard_interp.mat"
  "data_ipics/UNC_Martos_interp.mat"
  "data_ipics/UNC_Shen_interp.mat"
  "data_ipics/UNC_Stal_interp.mat"
  "data_ipics/BMIN_Losing_interp.mat"
  "data_ipics/BMAX_Losing_interp.mat"
  "data_ipics/specularity.mat"
  "data_ipics/coldex_icethk.mat"
  "Ts_interp.mat"
  "Mb_interp.mat"
  "data_ipics/mouginot_icevel.mat"
  "data_ipics/sink_mask_new.mat"
)

# loop over files and copy if they exist
for f in "${FILES[@]}"; do
    if [ -f "$f" ]; then
        cp -v "$f" "$DEST/"
    else
        echo "WARNING: File not found: $f"
    fi
done

echo "✅ Copy complete. Files saved in $DEST/"

