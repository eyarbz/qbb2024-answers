#!/bin/bash

echo -e "MAF\tFeature\tEnrichment" > snp_counts.txt

for MAF_value in 0.1 0.2 0.3 0.4 0.5
do
    MAF_file=chr1_snps_${MAF_value}.bed
    bedtools coverage -a genome_chr1.bed -b ${MAF_file} > chr1_coverage_${MAF_value}.txt
    snp_sum_genome=$(awk '{s+=$4}END{print s}' chr1_coverage_${MAF_value}.txt)
    chr_size=$(awk '{s+=$6}END{print s}' chr1_coverage_${MAF_value}.txt)
    background=$(bc -l -e "${snp_sum_genome}/${chr_size}")
    for feature in genes exons cCREs introns other
    do
        feature_file=${feature}_chr1.bed
        bedtools coverage -a ${feature_file} -b ${MAF_file} > ${feature_file}_coverage_${MAF_value}.txt
        snp_sum_feature=$(awk '{s+=$4}END{print s}' ${feature_file}_coverage_${MAF_value}.txt)
        feature_size=$(awk '{s+=$6}END{print s}' ${feature_file}_coverage_${MAF_value}.txt)
        ratio=$(bc -l -e "${snp_sum_feature}/${feature_size}")
        enrichment=$(bc -l -e "${ratio}/${background}")
        echo -e "${MAF_value}\t${feature}\t${enrichment}" >> snp_counts.txt
    done
done


    
