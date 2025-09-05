#!/bin/bash

# -------------------------------
# Script to run a single simulation chunk
# Usage: ./sim_one_file.sh <chunk_name> <merged_tsv> <output_root> <log_file>
# -------------------------------

# Arguments
chunk_name=$1       # Name of the current chunk
merged_tsv=$2       # Path to the merged TSV input file
output_root=$3      # Path to the output ROOT file
log_file=$4         # Path to the log file

# Log the start of the simulation chunk
echo "[$(date)] --- Starting chunk $chunk_name ---" | tee -a "$log_file"
echo "🔧 Running simulation for $chunk_name"

# Load Geant4 environment
source /home/dominik/G4Simulation/profileG4.sh

# Run the custom simulation
/home/dominik/G4Simulation/build/cmd/custom_simulation "$output_root" \
    -det 226:55:2 \
    -mask 467:170:128.3:101.4:20 \
    -n 1 \
    -masktype nowallpetcut \
    -cutx 57 \
    -cuty 45 \
    -sourceBins 140:200 \
    -nlay 7 \
    -source 0:0 \
    -1d \
    -sFile "$merged_tsv" \
    -fullscat

# Capture the exit code
exit_code=$?

# Check if the simulation succeeded
if [ $exit_code -ne 0 ]; then
    # Check if the process was killed by a signal (>128)
    if [ $exit_code -gt 128 ]; then
        signal=$((exit_code - 128))
        echo "⚠️ Simulation $chunk_name crashed with signal $signal at $(date)" | tee -a "$log_file"
    else
        echo "⚠️ Simulation $chunk_name failed with exit code $exit_code at $(date)" | tee -a "$log_file"
    fi
    echo "Skipping this chunk." | tee -a "$log_file"
else
    echo "✅ Simulation $chunk_name completed successfully at $(date)" | tee -a "$log_file"
    echo "Output ROOT file: $output_root" | tee -a "$log_file"
fi
