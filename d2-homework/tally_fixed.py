#!/usr/bin/env python3

import sys

my_file = open( sys.argv[1] )

chr = ""
count = 0

for my_line in my_file:

    if "#" in my_line:
        continue
        
    fields = my_line.split("\t")

    if fields[0] != chr:

        if chr != "":
            print( count, chr )

        chr = fields[0]
        count = 1 #changed from 0
        continue

    else:
        count = count + 1

    #count = count + 1

print(count, chr) #added this to 

my_file.close()