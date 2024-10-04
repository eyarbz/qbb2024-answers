#!/usr/bin/env python3

with open('biallelic.vcf', 'r') as infile:
    allele_freqs = []
    read_depths = []
    for line in infile:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        #allele frequency spectrum
        info = fields[7] #info is in 7th field (0 based)
        AF = info.split(';')[3] #allele frequency (AF) is in 3rd semicolon delimited field (0 based)
        AF = AF[3:] #remove 'AF=' from string
        allele_freqs.append(AF)
        #read depth distribution
        format = fields[9] #format is in 9th field (0 based)
        DP = format.split(':')[2]
        read_depths.append(DP)

with open('AF.txt', 'w') as outfile:
    outfile.write('Allele_Frequency'+'\n')
    for AF in allele_freqs:
        outfile.write(f"{AF}\n")

with open('DP.txt', 'w') as outfile:
    outfile.write('Read_Depth'+'\n')
    for DP in read_depths:
        outfile.write(f"{DP}\n")