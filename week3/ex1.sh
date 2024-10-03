#!/usr/bin/env bash

### Question 1.1 ###
read=$(sed -n '2p' A01_09.fastq)
length=${#read}
echo ${length}
#76 bases

### Question 1.2 ###
lines_tot=$(wc -l < A01_09.fastq)
reads_tot=$(bc -l -e "${lines_tot}/4")
echo ${reads_tot}
#669548 reads

### Question 1.3 ###
genome_len=12100000
coverage=$(bc -l -e "(${length}*${reads_tot})/${genome_len}")
echo ${coverage}
#Coverage: 4.2x

### Question 1.4 ###
fastq_sizes=$(du -h *.fastq)
echo ${fastq_sizes}
#Max: 149M A01_62.fastq
#Min: 110M A01_27.fastq

### Question 1.5 ###
#fastqc *.fastq
#The median looks like it is about 36
#Probability is 10^(-Q(A)/10), subbing Q(A) = 36 gives a low probability of 0.00025
#Not too much variation, boxes are consistent size throughout length of read. Outliers do seem to get less extreme to the center of the read.