#!/usr/bin/env python3

import sys


file = open(sys.argv[2])

for line in file:
    line = line.rstrip("\n")
    if sys.argv[1] in line:
        print(line)

file.close()