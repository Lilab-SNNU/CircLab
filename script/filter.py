#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/3/8 下午6:24
# @Author  : xm@10.150.62.165
# @File    : filter.py
# @Description :
import pandas as pd
import sys

def main():
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    threshold = snakemake.config['threshold']

    df = pd.read_csv(input_file, sep='\t')
    df['count'] = df['Source'].apply(lambda x: x.count(','))

    if threshold == 1:
        filtered_df = df
    elif threshold == 2:
        filtered_df = df[df['count'] >= 1]
    elif threshold == 3:
        filtered_df = df[df['count'] >= 2]
    elif threshold == 4:
        filtered_df = df[df['count'] >= 3]
    else:
        raise ValueError(f"Invalid threshold: {threshold} (must be 1-4)")

    filtered_df.drop('count', axis=1).to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()