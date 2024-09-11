#!/usr/bin/env python3

import sys

import numpy as np

#1: gene_tissue_tsv
#Output: Dictionary where gene ID is key and tissue is value

file_tsv = sys.argv[1]

with open(file_tsv, mode = 'r') as fs_tsv:
    geneid_tissue_dict = {}

    for line in fs_tsv:
        fields = line.strip('\n').split('\t')
        geneid_tissue_dict[fields[0]] = fields[2]


#print(geneid_tissue_dict)

#2: GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt (metadata)
#Output: Dictionary where key is tissue and value is list of sample IDs corresponding to that tissue

file_metadata = sys.argv[2]

with open(file_metadata, mode = 'r') as fs_metadata:
    tissue_sampid_dict = {}

    fs_metadata.readline()

    for line in fs_metadata:
        fields = line.strip('\n').split('\t')
        #make the key the tissue (column 6 of metadata file)
        key = fields[6]
        #make the value the sample id (column 0 of metadata file)
        value = fields[0]
        #for each new tissue, make the value a list
        tissue_sampid_dict.setdefault(key, [])
        #append the sample id to the value list
        tissue_sampid_dict[key].append(value)


#print(tissue_sampid_dict)

#3: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct (gene expression file)
#Output: List of sample IDs

file_gene_exp = sys.argv[3]

with open(file_gene_exp, mode = 'r') as fs_gene_exp:
    #read the through the first 2 lines of the file
    fs_gene_exp.readline()
    fs_gene_exp.readline()
    #make a list of the samp ids (column headers starting from column 2)
    samp_ids = fs_gene_exp.readline().strip('\n').split('\t')[2:]

#print(samp_ids)

#4: Creating dictionary of sample ids and their column indices from gene expression file
#Output: (Relating to 3) Dictionary where key is sample ID and value is index of sample ID in the list

samp_ids_index_dict = {}

for samp_id in samp_ids:
    samp_ids_index_dict[samp_id] = samp_ids.index(samp_id)

#print(samp_ids_index_dict)

#5: Create the your dictionary of each tissue and which column indices in the gene expression file correspond to it

tissue_index_dict = {}

#loop through each tissue from the tissue_sampid_dict (key = tissue, value = list of sample ids)
for key in tissue_sampid_dict:
    #make each dictionary value a list
    tissue_index_dict.setdefault(key, [])

    #for each list value in tissue_sampid_dict, loop through the sample ids in that list
    for samp_id in tissue_sampid_dict[key]:
        
        #find which column index corresponds to the sample id loop variable and append that index to the list value in tissue_index_dict
        if samp_id in samp_ids_index_dict:
            tissue_index_dict[key].append(samp_ids_index_dict[samp_id] )
            
        
#print(tissue_index_dict)

#Figure out which tissues have the highest and lowest number of samples

#Create a dictionary to hold number of samples as key and tissue as value
no_sample_tissue_dict = {}

for key in tissue_index_dict:
    value = tissue_index_dict[key]
    no_samples = len(value)
    no_sample_tissue_dict[no_samples] = key

keys = no_sample_tissue_dict.keys()
max_samples = max(keys)
min_samples = min(keys)

#print(no_sample_tissue_dict)
#print("Maximum samples:", no_sample_tissue_dict[max_samples], max_samples)
#print("Minimum samples:", no_sample_tissue_dict[min_samples], min_samples)

#Max: Muscle - Skeletal 803 samples
#Min: Cells - Leukemia cell line (CML) 0


#6: Steps:
#a: From expression file, look at each gene id (rows) and figure out if it is in the gene id: tissue dictionary.
#b: If it is, record it, the corresponding tissue, and figure out which columns in the expression file correspond to the tissue.
#b_1: To do this, find the list of indices that correspond to said tissue from the tissue_index_list dictionary (key = tissue, value = list of indices)
#b_2: Then, index the row using that list of indices to pull out expressions from relevant tissues
#c: Add gene id, tissue, and expression values to their own lists

file_gene_exp = sys.argv[3]

with open(file_gene_exp, mode = 'r') as fs_gene_exp:
    #read the through the first 2 lines of the file
    fs_gene_exp.readline()
    fs_gene_exp.readline()
    fs_gene_exp.readline()

    gene_id_list = []
    tissue_list = []
    expressions_list = []

    for line in fs_gene_exp:
        fields = line.strip('\n').split('\t')
        gene_id = fields[0]
        if gene_id in geneid_tissue_dict:
            tissue = geneid_tissue_dict[gene_id]
            indices = tissue_index_dict[tissue]
            expressions = np.array(fields[2:])
            tissue_expressions = expressions[indices]

            gene_id_list.append(gene_id)
            tissue_list.append(tissue)
            expressions_list.append(tissue_expressions)

#print(gene_id_list, tissue_list, expressions_list)


#7: Creating tsv file with gene id, tissue, and expression values

tsv_path = '/Users/cmdb/qbb2024-answers/d4-lunch/id_tissue_expressions.tsv'

with open(tsv_path, 'w') as fh:
    headers = "gene_id" + "\t" + "tissue" + "\t" + "expression" + "\n"
    fh.write(headers)

    for i in range(len(gene_id_list)):
        geneid = gene_id_list[i]
        tissue = tissue_list[i]
        expression_array = expressions_list[i]

        for j in range(len(expression_array)):
            expression = expression_array[j]
            line = f"{geneid}\t{tissue}\t{expression}\n"
            fh.write(line)