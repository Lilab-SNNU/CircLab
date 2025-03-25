#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
@Project ：_circlab 
@File    ：conservation_analysis.py.py
@Author  ：xm
@Date    ：2024/12/30 下午3:50 
'''
import os
import argparse
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import combinations
from Bio.SeqIO import parse
import pandas as pd

def parse_species_list(species_file):
    with open(species_file, 'r') as f:
        return [line.strip() for line in f]

def load_circRNA_sequences(folder, species_list):
    circRNA_data = {}
    for species in species_list:
        file_path = os.path.join(folder, f"{species}.fa")
        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} not found, skipping.")
            continue

        sequences = {}
        for record in parse(file_path, "fasta"):
            sequences[record.id] = str(record.seq)
        circRNA_data[species] = sequences
    return circRNA_data

def pairwise_compare(seq_dict1, seq_dict2, aligner, score_threshold=50, flank=25):
    conserved = set()
    for circ1, seq1 in seq_dict1.items():
        extended_seq1 = seq1[:flank] + seq1 + seq1[-flank:]
        for circ2, seq2 in seq_dict2.items():
            extended_seq2 = seq2[:flank] + seq2 + seq2[-flank:]

            seq_record1 = SeqRecord(Seq(extended_seq1), id=circ1)
            seq_record2 = SeqRecord(Seq(extended_seq2), id=circ2)

            alignments = aligner.align(seq_record1.seq, seq_record2.seq)
            best_alignment = max(alignments, key=lambda x: x.score)

            if best_alignment.score >= score_threshold:
                conserved.add(circ1)
    return conserved

def build_conservation_matrix(circRNA_data, species_list, aligner, score_threshold=50, flank=25):
    all_circRNAs = set()
    pairwise_results = {}

    for species1, species2 in combinations(species_list, 2):
        seq_dict1 = circRNA_data.get(species1, {})
        seq_dict2 = circRNA_data.get(species2, {})
        conserved = pairwise_compare(seq_dict1, seq_dict2, aligner, score_threshold, flank)
        pairwise_results[(species1, species2)] = conserved
        all_circRNAs.update(conserved)

    matrix = {circ: {species: 0 for species in species_list} for circ in all_circRNAs}
    for (species1, species2), conserved_set in pairwise_results.items():
        for circ in conserved_set:
            matrix[circ][species1] = 1
            matrix[circ][species2] = 1
    return matrix

def find_fully_conserved_circRNAs(matrix, species_list):

    fully_conserved = [
        circ for circ, species_status in matrix.items()
        if all(species_status[species] == 1 for species in species_list)
    ]
    return fully_conserved

def main(species_file, input_folder, output_file, score_threshold=50, flank=25):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aligner.match_score = 1
    aligner.mismatch_score = -1

    species_list = parse_species_list(species_file)

    circRNA_data = load_circRNA_sequences(input_folder, species_list)

    conservation_matrix = build_conservation_matrix(
        circRNA_data, species_list, aligner, score_threshold, flank
    )

    fully_conserved = find_fully_conserved_circRNAs(conservation_matrix, species_list)

    with open(output_file, 'w') as f:
        for circ in fully_conserved:
            f.write(f"{circ}\n")
    print(f"Fully conserved circRNAs written to {output_file}")

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="Cross-species conservation analysis of circRNAs")
    parser.add_argument("-s", "--species", required=True, help="File containing species names (species.txt)")
    parser.add_argument("-i", "--input_folder", required=True, help="Folder containing circRNA sequences for each species")
    parser.add_argument("-o", "--output_file", required=True, help="Output file for conserved circRNAs")
    parser.add_argument("-t", "--score_threshold", type=int, default=50, help="Alignment score threshold (default: 50)")
    parser.add_argument("-f", "--flank", type=int, default=25, help="Length of upstream and downstream flanking sequences (default: 25bp)")

    args = parser.parse_args()
    main(args.species, args.input_folder, args.output_file, args.score_threshold, args.flank)
