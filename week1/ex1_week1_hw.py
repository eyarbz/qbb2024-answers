#!/usr/bin/env python3

import numpy as np
import scipy as sp

##Exercise 1

###Simulation:
genome_size = 1000000
coverage = 30
read_size = 100
output_path = "/Users/cmdb/qbb2024-answers/week1/genome_coverage.txt"
#Generate total number of reads based on genome size, desired coverage, and read size
num_reads = int((genome_size * coverage)/read_size)
#Initialize a genome with no mapped reads to any bases
genome_coverage = np.zeros(genome_size)

for _ in range(num_reads):
    #Generate a random base in the genome for a mapped read to begin at
    start_pos = np.random.randint(0, genome_size - read_size + 1)
    #Generate the end of that mapped read
    end_pos = start_pos + read_size
    #Add a count of 1 to the region a read maps to in the genome
    genome_coverage[start_pos:end_pos] += 1

#Save genome coverage to output txt file
np.savetxt(output_path, genome_coverage, fmt='%d', header='Coverage', comments='')


###Poisson Distribution:
output_path_poisson = "/Users/cmdb/qbb2024-answers/week1/poisson.txt"
#Find the highest coverage from the simulation
max_coverage = int(max(genome_coverage))
#Generate a range of coverages from the simulation
read_list = np.arange(0, max_coverage + 1)
#Initialize an array to store the frequencies of finding each coverage from poisson distribution
poisson_freq_list = np.zeros(len(read_list))

for i in read_list:
    poisson_freq_list[i] = sp.stats.poisson.pmf(i , mu = coverage) * genome_size

#Put coverage number in one column, and the frequencies in the next
poisson_df = np.column_stack((read_list, poisson_freq_list))
#Save poisson distributed coverages to a txt file
np.savetxt(output_path_poisson, poisson_df, fmt='%d', header='Coverage Frequency', comments='')


###Normal Distribution
output_path_normal = "/Users/cmdb/qbb2024-answers/week1/normal.txt"
std = np.sqrt(coverage)
#Initialize an array to store the frequencies of finding each coverage from normal distribution
normal_freq_list = np.zeros(len(read_list))

for i in read_list:
    normal_freq_list[i] = sp.stats.norm.pdf(i, coverage, std) * genome_size

#Put coverage number in one column, and the frequencies in the next
normal_df = np.column_stack((read_list, normal_freq_list))
#Save normally distributed coverages to a txt file
np.savetxt(output_path_normal, normal_df, fmt='%d', header='Coverage Frequency', comments='')