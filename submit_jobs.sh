#!/bin/bash
h_values=()
for h in $(seq 0.0 0.1 5); do
   formatted_h=$(printf "%.2f" "$h")
   h_values+=("$formatted_h")
done

# Loop over each h value
for h in "${h_values[@]}"; do
    # Create a unique job script for each h value
    job_file="SMM_h_${h}.run"
    sed "s/\${h_val}/$h/g" SSM.run > "$job_file"

    # Submit the job
    sbatch -A ctmc "$job_file"

    echo "Submitted job for h=$h"
done