#!/bin/bash

# -------------------------------
# Script to run simulations on TSV files in chunks
# Converts chunks to ROOT files using sim_one_file.sh
# Supports parallel execution with limited concurrent jobs
# -------------------------------

# Directory containing filtered TSV files
files_dir="/scratch3/gccb/simulations/PhaseSpaceFiles/additional_energies/extracted_tsv_filtered"

# Output directory for simulation results
sim_dir="/scratch3/gccb/simulations/PhaseSpaceFiles/additional_energies/sim_fixed_geo"

# Number of TSV files to process in one chunk
chunk_size=10

# Log file for the simulation process
log_file="${sim_dir}/simulation_log.txt"

# Create output directory if it does not exist
mkdir -p "$sim_dir"

# Enable nullglob to prevent errors if no files match a pattern
shopt -s nullglob

# Maximum number of parallel jobs
max_jobs=6
job_count=0

# Optional: Function to release memory cache (requires sudo)
release_cache() {
    echo "🧹 Releasing cache..."
    sudo sync
    sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
}

# Log start time
echo "=== Simulation started at $(date) ===" >> "$log_file"

# Loop over all subdirectories in the TSV files folder
for dir in "$files_dir"/*/; do
    dir_name=$(basename "$dir")
    output_subdir="${sim_dir}/${dir_name}"
    mkdir -p "$output_subdir"

    echo "📂 Processing directory: $dir_name"

    # Get all TSV files in this directory
    tsv_files=("$dir"/*.tsv)
    total_files=${#tsv_files[@]}

    # Loop over TSV files in chunks of size $chunk_size
    for (( i=0; i<total_files; i+=chunk_size )); do
        (
            # Extract current chunk
            chunk_files=("${tsv_files[@]:i:chunk_size}")
            chunk_name="${dir_name}_chunk_$((i/chunk_size + 1))"
            merged_tsv="${output_subdir}/${chunk_name}.tsv"
            output_root="${output_subdir}/${chunk_name}.root"

            # Skip chunk if merged TSV already exists
            if [[ -f "$merged_tsv" ]]; then
                echo "⚠️ Chunk $chunk_name already exists. Skipping simulation." | tee -a "$log_file"
                exit 0
            fi

            # Log chunk start
            echo "[$(date)] --- Starting chunk $chunk_name ---" | tee -a "$log_file"
            echo "Files in chunk:" >> "$log_file"
            for f in "${chunk_files[@]}"; do
                echo "  $f" >> "$log_file"
            done

            # Merge TSV files in the chunk
            echo "📦 Merging ${#chunk_files[@]} files into $merged_tsv"
            # Skip header lines except for the first file
            awk 'FNR==1 && NR!=1 { next } { print }' "${chunk_files[@]}" > "$merged_tsv"

            # Run simulation script on merged TSV
            ./sim_one_file.sh "$chunk_name" "$merged_tsv" "$output_root" "$log_file"
        ) &  # Run chunk in background

        # Manage parallel jobs
        ((job_count++))
        if (( job_count % max_jobs == 0 )); then
            wait  # Wait for all background jobs to finish
        fi
    done
done

# Wait for any remaining background jobs to finish
wait

# Log completion time
echo "=== All simulations completed at $(date) ===" | tee -a "$log_file"
