#!/usr/bin/env python3

import sys

file = open(sys.argv[1])

for line in file:
    if "##" in line:
        continue
    
    fields = line.split("\t")
    chrom = fields[0]
    start = fields[3]
    stop = fields[4]

    thing = fields[-1]
    fields2 = thing.split(";")
    
    gene = ""

    for entry in fields2:
        if "gene_name" in entry:
            gene = entry
        
    gene = gene.rstrip("\"")
    gene = gene.lstrip("gene_name \"")

    print(chrom + "\t" + start + "\t" + stop + "\t" + gene)



file.close()