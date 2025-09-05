#!/bin/bash

# -------------------------------
# Script to merge ROOT files in subdirectories using hadd
# Skips small files, limits combined file size to max_chunk_size
# -------------------------------

# Input directory containing simulation ROOT files
input_dir="/scratch3/gccb/simulations/PhaseSpaceFiles/additional_energies/sim_fixed_geo"

# Output directory for combined ROOT files
output_dir="${input_dir}/combined"

# Log file for process output
log_file="${output_dir}/combine_log.txt"

# Minimum file size to consider (in bytes)
min_size=50000   # 50 KB

# Maximum total size for a single combined file (in bytes)
max_chunk_size=$((300 * 1024 * 1024))  # 300 MB

# Create output directory if it does not exist
mkdir -p "$output_dir"

# Initialize log file
echo "📝 Starting ROOT files merge..." > "$log_file"

# Loop over all subdirectories in the input directory
for subdir in "$input_dir"/*; do
    # Skip the combined output directory itself
    if [[ "$(basename "$subdir")" == "combined" ]]; then
        continue
    fi

    sub_name=$(basename "$subdir")
    echo "📁 Processing $sub_name" | tee -a "$log_file"

    # Output path for the combined ROOT file
    combined_output="${output_dir}/${sub_name}_combined.root"

    # Array to hold valid ROOT files to merge
    root_files=()
    current_size=0  # Track cumulative size of files to merge

    # Loop over all ROOT files in the current subdirectory
    for root_file in "$subdir"/*.root; do
        # Skip if no files match
        if [ ! -f "$root_file" ]; then continue; fi

        # Get file size in bytes
        size=$(stat -c%s "$root_file")

        # Skip files smaller than minimum threshold
        if [ "$size" -lt "$min_size" ]; then
            echo "❌ Skipped (too small): $root_file ($size B)" | tee -a "$log_file"
            continue
        fi

        # Stop adding files if max chunk size would be exceeded
        if (( current_size + size > max_chunk_size )); then
            echo "⏹️ Stopping: max chunk size reached." | tee -a "$log_file"
            break
        fi

        # Add file to array and update cumulative size
        root_files+=("$root_file")
        current_size=$((current_size + size))

        echo "✅ Added: $root_file ($size B)" | tee -a "$log_file"
    done  # end inner loop over ROOT files

    # Merge valid ROOT files if any
    if [ "${#root_files[@]}" -gt 0 ]; then
        echo "🔧 Running hadd for $sub_name..." | tee -a "$log_file"
        hadd -f "$combined_output" "${root_files[@]}" >> "$log_file" 2>&1
        echo "🎯 Combined output: $combined_output" | tee -a "$log_file"
    else
        echo "⚠️ No valid files to combine in $sub_name" | tee -a "$log_file"
    fi

    echo "--------------------------------------------------" >> "$log_file"
done  # end outer loop over subdirectories

echo "✅ All combinations complete. Log saved to: $log_file"
