#!/usr/bin/env bash

bwa index sacCer3.fa

for my_sample in *.fastq
do
    my_sample=`basename ${my_sample} .fastq`
    bwa mem -R "@RG\tID:${my_sample}\tSM:${my_sample}" sacCer3.fa ${my_sample}.fastq > ${my_sample}.sam
    samtools sort -@ 4 -O bam -o ${my_sample}.bam ${my_sample}.sam
    samtools index ${my_sample}.bam
done

### Question 2.2 ###
tail -n +21 A01_09.sam | wc -l
#There are 669548 reads

### Question 2.3 ###
tail -n +21 A01_09.sam | cut -d $'\t' -f 3 | grep "chrIII" | wc -l
#17815 alignments are to loci on chromosome 3

### Question 2.4 ###
#It is a little hard to tell, it looks like the coverage is a little higher, but there is a lot of variation

### Question 2.5 ###
#I see 3 SNPs. I am uncertain about the rightmost SNP (chrI:113,326): only 2 reads show the SNP, and one SNP is at pos1 of read, which may not be great quality.

### Question 2.6 ###
#SNP is at chrIV:825,834. It does not fall within a gene.