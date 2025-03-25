#!/bin/bash
# Parameters passed by Snakemake
bwa_name_file="$1"
bwa_result_file="$2"
star_name_file="$3"
star_result_file="$4"
ciri_name_file="$5"
ciri_result_file="$6"
findcirc_name_file="$7"
findcirc_result_file="$8"

# Function to process BWA and STAR files (same format)
process_bwa_star_files() {
    local name_file=$1
    local result_file=$2
    local circ_files=("$3")

    > "$result_file"

    for file in "${circ_files[@]}"; do
        while IFS=$'\t' read -r chrom start end; do
            awk -v chrom="$chrom" -v start="$start" -v end="$end" \
                'BEGIN {FS=OFS="\t"} $1 == chrom && $2 == start && $3 == end {print $0}' "$file" >> "$result_file"
        done < "$name_file"
    done

    awk '!seen[$0]++' "$result_file" > temp && mv temp "$result_file"
}

# Function to process CIRI files (different format)
process_ciri_files() {
    local name_file=$1
    local result_file=$2
    local circ_files=("$3")

    > "$result_file"

    for file in "${circ_files[@]}"; do
        while IFS=$'\t' read -r chrom start end; do
            awk -v chrom="$chrom" -v start="$start" -v end="$end" \
                'BEGIN {FS=OFS="\t"} $2 == chrom && $3 == start && $4 == end {print $0}' "$file" >> "$result_file"
        done < "$name_file"
    done

    awk '!seen[$0]++' "$result_file" > temp && mv temp "$result_file"
}

# Assuming your script gets arrays of files if multiple samples per type may exist:
bwa_circ_files=( $(find /path/to/bwa/files -name 'circ*.txt') )
star_circ_files=( $(find /path/to/star/files -name 'circ*.txt') )
ciri_circ_files=( $(find /path/to/ciri/files -name 'ciri_*.txt') )
# Placeholder: findcirc_circ_files also to be specified.

# Process BWA and STAR files
process_bwa_star_files "$bwa_name_file" "$bwa_result_file" "${bwa_circ_files[@]}"
process_bwa_star_files "$star_name_file" "$star_result_file" "${star_circ_files[@]}"

# Process CIRI files
process_ciri_files "$ciri_name_file" "$ciri_result_file" "${ciri_circ_files[@]}"

# If FindCirc also needs processing, add a similar function call
# process_findcirc_files "$findcirc_name_file" "$findcirc_result_file" "${findcirc_circ_files[@]}"