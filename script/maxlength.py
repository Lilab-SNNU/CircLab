#!/usr/bin/env python
'''
Modify according to find_circ maxlength.py
'''
import sys

max_l = float(sys.argv[1])

for line in sys.stdin:
    if line.startswith("#"):
        print(line.rstrip())
        continue
    
    cols = line.rstrip().split("\t")
    try:
        start, end = float(cols[1]), float(cols[2])
        if end - start <= max_l:
            print(line.rstrip())
    except ValueError:
         continue 
