#!/usr/bin/env python3

##Exercise 2:

reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']
k = 3
graph = set()
output_path = '/Users/cmdb/qbb2024-answers/week1/ex2_edges.dot'

for read in reads:
    for i in range(len(read) - k):
        kmer1 = read[i: i+k]
        kmer2 = read[i+1: i+1+k]
        edge = f"{kmer1} -> {kmer2}"
        graph.add(edge)
    
with open(output_path, 'w') as fh:
    fh.write('digraph {' + '\n')
    for edge in graph:
        fh.write(edge + '\n')
    fh.write('}')

print(f"There are {len(graph)} edges")