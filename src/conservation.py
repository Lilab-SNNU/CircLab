#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
@Project ：_circlab
@File    ：conservation.py
@Author  ：xm
@Date    ：2024/12/30 下午3:50
'''

import sys
import os
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import combinations
from Bio.SeqIO import parse
import pandas as pd
import argparse
from multiprocessing import Pool

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
            new_id = f"{record.id}|{species}"
            sequences[new_id] = str(record.seq)
        circRNA_data[species] = sequences
    return circRNA_data

def pairwise_compare(args):
    species1, seq_dict1, species2, seq_dict2, aligner, score_threshold, flank = args
    conserved = set()
    conserved_details = {}

    required_length = 2 * flank
    total_pairs = len(seq_dict1) * len(seq_dict2)
    processed_pairs = 0

    for circ1, seq1 in seq_dict1.items():
        if len(seq1) < required_length:
            continue
        extended_seq1 = (seq1[-flank:] + seq1[:flank]).upper()

        for circ2, seq2 in seq_dict2.items():
            if len(seq2) < required_length:
                continue  

            extended_seq2 = (seq2[-flank:] + seq2[:flank]).upper()

            seq_record1 = SeqRecord(Seq(extended_seq1), id=circ1)
            seq_record2 = SeqRecord(Seq(extended_seq2), id=circ2)

            alignments = aligner.align(seq_record1.seq, seq_record2.seq)
            best_alignment = max(alignments, key=lambda x: x.score)

            if best_alignment.score >= score_threshold:
                conserved.add(circ1)
                if circ1 not in conserved_details:
                    conserved_details[circ1] = set()
                conserved_details[circ1].add(circ2)

            processed_pairs += 1
            if processed_pairs % 5000 == 0 or processed_pairs == total_pairs:
                print(f"    Compared {species1} vs {species2}: {processed_pairs}/{total_pairs} pairs processed...")

    print(f"  Compared {species1} vs {species2}: Found {len(conserved)} conserved circRNAs.")
    return (species1, species2, conserved, conserved_details)


def build_conservation_matrix(circRNA_data, species_list, aligner, score_threshold=50, flank=25, processes=4):
    all_circRNAs = set()
    pairwise_results = {}

    args_list = []
    for species1, species2 in combinations(species_list, 2):
        seq_dict1 = circRNA_data.get(species1, {})
        seq_dict2 = circRNA_data.get(species2, {})
        args_list.append((species1, seq_dict1, species2, seq_dict2, aligner, score_threshold, flank))

    with Pool(processes=processes) as pool:
        results = pool.map(pairwise_compare, args_list)

    for species1, species2, conserved_set, conserved_details in results:
        pairwise_results[(species1, species2)] = conserved_details
        all_circRNAs.update(conserved_set)

    matrix = {circ: {species: 0 for species in species_list} for circ in all_circRNAs}
    conserved_names_mapping = {circ: set() for circ in all_circRNAs}

    for (species1, species2), conserved_details in pairwise_results.items():
        for circ1, conserved_circs in conserved_details.items():
            matrix[circ1][species1] = 1
            matrix[circ1][species2] = 1
            conserved_names_mapping[circ1].update(conserved_circs)

    return matrix, conserved_names_mapping


def find_fully_conserved_circRNAs(matrix, conserved_names_mapping, circRNA_data, species_list):
    conservation_stats = []
    for circ, species_status in matrix.items():
        try:
            source_species = circ.split('|')[-1]
        except IndexError:
            print(f"Invalid circRNA ID format: {circ}")
            continue

        conservation_count = sum(
            1 for species in species_list if species != source_species and species_status[species] > 0
        )

        conserved_circRNAs = [
            f"{circ2}" for circ2 in conserved_names_mapping[circ]
            for species in species_list if circ2 in circRNA_data[species]
        ]

        conserved_circRNA_names = ";".join(conserved_circRNAs)

        circ_seq = next(
            (circRNA_data[species][circ] for species in species_list if species != source_species and circ in circRNA_data[species]),
            ""
        )

        if not circ_seq:
            circ_seq = next(
                (circRNA_data[species][circ] for species in species_list if circ in circRNA_data[species]),
                ""
            )

        conservation_stats.append((circ, conservation_count, circ_seq, conserved_circRNA_names))

    conservation_stats.sort(key=lambda x: x[1], reverse=True)
    return conservation_stats


def save_conservation_matrix(matrix, output_file):
    df = pd.DataFrame.from_dict(matrix, orient='index')
    df.to_csv(output_file, sep='\t')

def save_to_file(results, output_file):
    with open(output_file, 'w') as f:
        for circ, count, seq, conserved_names in results:
            f.write(f"{circ}\t{count}\t{seq}\t{conserved_names}\n")

def main(argv=None):
    parser = argparse.ArgumentParser(description="Conservation Analysis Module")

    parser.add_argument("-s", "--species", required=True, help="File containing species names (species.txt)")
    parser.add_argument("-i", "--input_folder", required=True, help="Folder containing circRNA sequences for each species")
    parser.add_argument("-o", "--output_file", required=True, help="Output file for conservation matrix")
    parser.add_argument("-t", "--score_threshold", type=int, default=50, help="Alignment score threshold (default: 50)")
    parser.add_argument("-f", "--flank", type=int, default=25, help="Length of flanking regions to include (default: 25bp)")
    parser.add_argument("-p", "--processes", type=int, default=4, help="Number of processes for parallel computation (default: 4)")

    args = parser.parse_args(argv)

    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aligner.match_score = 1
    aligner.mismatch_score = -1

    species_list = parse_species_list(args.species)

    circRNA_data = load_circRNA_sequences(args.input_folder, species_list)

    conservation_matrix, conserved_names_mapping = build_conservation_matrix(
        circRNA_data, species_list, aligner, args.score_threshold, args.flank, args.processes
    )

    save_conservation_matrix(conservation_matrix, args.output_file)
    print(f"Conservation matrix written to {args.output_file}")

    conservation_stats = find_fully_conserved_circRNAs(conservation_matrix, conserved_names_mapping, circRNA_data, species_list)

    fully_conserved_file = args.output_file.replace('.txt', '_fully_conserved.txt')
    save_to_file(conservation_stats, fully_conserved_file)
    print(f"Fully conserved circRNAs written to {fully_conserved_file}")

if __name__ == "__main__":
    main(sys.argv[1:])



