#!/bin/bash
prefix="prefix_of_the_array_data" # prefix of the array data
rs_number="rs761167763"  #rsnumber of the variant

# the name of flanking region or snp
ranges=("100kb" "1000kb" "10000kb" "50kb" "20snp")

for range in "${ranges[@]}"; do
    echo "Processing range: $range"

    java -Xmx4g -jar /your/path/to/beagle.jar \
        gt=${prefix}_${rs_number}_${range}.vcf \
        out=${prefix}_${rs_number}_${range}_5.phased \
        gp=true \
        impute=true \
        nthreads=4
done

