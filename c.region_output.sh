#!/bin/bash

plink=/opt/ohpc/Taiwania3/pkg/biology/PLINK/PLINK_v1.90/plink
rs_number="rs747512450"  
prefix="merged_nofamily_nonosex.qc.noindel"

regions=("20snp" "50kb" "100kb" "1000kb")

for region in "${regions[@]}"; do
    pos_file="filtered_${rs_number}_${region}_positions.txt"
    range_file="${rs_number}_${region}_range.txt"
    
    awk 'NR > 1 {print $1"\t"$2"\t"$2}' "$pos_file" > "$range_file"

    $plink --bfile "$prefix" \
           --extract range "$range_file" \
           --make-bed \
           --out "${prefix}_${rs_number}_${region}"

    $plink --bfile "$prefix" \
           --extract range "$range_file" \
           --recode vcf \
           --out "${prefix}_${rs_number}_${region}"
done

