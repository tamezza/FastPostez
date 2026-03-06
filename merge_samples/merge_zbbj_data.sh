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

merge_sample "multijets_sherpa" 364686-364694
#merge_sample "multijets_364702" 364702
#merge_sample "multijets_364703" 364703
#merge_sample "multijets_364704" 364704
#merge_sample "multijets_364705" 364705
#merge_sample "multijets_364706" 364706
#merge_sample "multijets_364707" 364707
#merge_sample "multijets_364708" 364708
#merge_sample "multijets_364709" 364709
#merge_sample "multijets_364710" 364710
#merge_sample "multijets_364711" 364711
#merge_sample "multijets_364712" 364712
#merge_sample "top" 410470 410471 410644 410645 410658 410659
#merge_sample "WqqJets" 700843
#merge_sample "ZbbJets" 700855
#merge_sample "ZqqJets" 700849


hadd "$output_dir"/data.root "$input_dir"/ntuples/ntuples_1[5-8].root
