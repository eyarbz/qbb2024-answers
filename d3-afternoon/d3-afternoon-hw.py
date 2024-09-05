#!/usr/bin/env python3

import sys

import numpy as np 

#Question 1
fs = open(sys.argv[1], 'r')

fs.readline()
fs.readline()

fields = fs.readline()
fields = fields.strip("\n").split("\t")

#Lists for different elements of GTEx dataset
tissues = fields[2:]
gene_ids = []
gene_names = []
expression = []

for line in fs:
    fields = line.strip("\n").split("\t")
    gene_ids.append(fields[0]) 
    gene_names.append(fields[1])
    expression.append(fields[2:])

fs.close()

#Question 2

tissues = np.array(tissues)
gene_ids = np.array(gene_ids)
gene_names = np.array(gene_names)
expression = np.array(expression, dtype = float)
#print(expression)
#Because expression data has decimals

#Question 4
means = np.mean(expression[0:10], axis=1)
#print(means)


#Question 5
means_whole = np.mean(expression)
median_whole = np.median(expression)
#print(means_whole, median_whole)
#Mean ~ 16.55, median ~ 0.027
#Most genes have a low expression levels, but there are high outliers that skew the mean


#Question 6

log2_transform = np.log2(expression + 1)
means_log2 = np.mean(log2_transform)
medians_log2 = np.median(log2_transform)

#print(means_log2, medians_log2)
#Mean ~ 1.1, median ~ 0.039
#Mean has decreased a lot, median has not changed as significantly


#Question 7

log2_exp_sorted = np.sort(log2_transform, axis = 1)
diff_array = log2_exp_sorted[:,-1] - log2_exp_sorted[:,-2]


#Question 8
diff_array_sig = np.sum(diff_array >= 10)
#print(diff_array_sig)