#!/bin/bash

input_dir=$1
output_dir=$2
mkdir -p "$output_dir"

merge_sample() {
    local sample_name=$1
    shift
    local raw_args=("$@")
    local dsids=()

    # Expand any ranges (e.g., 364418-364422)
    for arg in "${raw_args[@]}"; do
        if [[ $arg =~ ^[0-9]+-[0-9]+$ ]]; then
            IFS='-' read -r start end <<< "$arg"
            for ((i=start; i<=end; i++)); do
                dsids+=("$i")
            done
        else
            dsids+=("$arg")
        fi
    done

    echo "Merging sample: $sample_name with DSIDs: ${dsids[*]}"

    # Collect matching files
    local file_list=()
    for dsid in "${dsids[@]}"; do
        matches=("$input_dir"/ntuples/*"${dsid}"*.root)
        if [ -e "${matches[0]}" ]; then
            file_list+=("${matches[@]}")
        else
            echo "Warning: No files found for DSID $dsid"
        fi
    done

    if [ "${#file_list[@]}" -gt 0 ]; then
        hadd -f "$output_dir/${sample_name}.root" "${file_list[@]}"
    else
        echo "No input files found for sample: $sample_name. Skipping."
    fi
}

#merge_sample "Zqqy_MG" 364418-364422
#merge_sample "gamma_jets_sherpa" 
merge_sample "tty" 545023-545024
merge_sample "Zqqy" 701287
merge_sample "Zbby" 701286
merge_sample "Wqqy" 701288
merge_sample "gamma_jets_pythia" 801667-801676

hadd "$output_dir"/data.root "$input_dir"/ntuples/ntuples_2[2-4].root


