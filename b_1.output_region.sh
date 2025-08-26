#!/bin/bash
plink=/your/path/to/PLINK/PLINK_v1.90/plink 
rs_number="rs761167763"  #rsnumber of the variant
prefix="prefix_of_the_array_data" #prefix of the array data

# the name of flanking region or snp
regions=("20snp" "50kb" "100kb" "1000kb" "10000kb")

for region in "${regions[@]}"; do
    range_file="${rs_number}_${region}_list.txt"

    # BED output
    $plink --bfile "$prefix" \
           --extract "$range_file" \
           --make-bed \
           --out "${prefix}_${rs_number}_${region}"

    # VCF output
    $plink --bfile "$prefix" \
           --extract "$range_file" \
           --recode vcf \
           --out "${prefix}_${rs_number}_${region}"
done

