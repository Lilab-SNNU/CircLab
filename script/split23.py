#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
for split ""
@Project ：pythonT 
@File    ：split23.py
@Author  ：xm
@Date    ：2023/10/23 上午10:08 
'''
f = open("annotation0.txt","r")
fo = open("annotation.txt",'w')
result = ''
for line in f:
    line1 = line.strip().split("\t")
    line_start = "\t".join(line1[:9])
    line_end = "\t".join(line1[11:])
    #print(line_start)
    start_site = line1[9].strip('""')
    end_site = line1[10].strip('""')
    #print(start_site)
    result += f'{line_start}\t{start_site}\t{end_site}\t{line_end}\n'
fo.write(result)
