#!/usr/bin/env python

"""
Script to generate CSV file containing each unique filepair from among provided FASTA files in order to blast them
using DIAMOND within Nextflow workflow.
"""

import argparse
import itertools
import csv
import sys
import os

from utils import get_full_filepaths, remove_files, get_filename, check_output_path


def parse_args(args=None):
    Description = "Create CSV file listing each unique pair of FASTA files to BLAST in (genome1, genome2) format."
    Epilog = "Example usage: python make_genome_filepairs.py <FASTA_PATH> <BLAST_PATH> <OUTPUT_PATH>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FASTA_PATH', metavar='fasta_path', type=str,
                        help='Path to directory containing neighborhood FASTA files.')
    parser.add_argument('BLAST_PATH', metavar='blast_path', type=str,
                        help='Path to directory containing BLAST text files of BLASTed neighborhoods.')
    parser.add_argument('OUTPUT_PATH', metavar='output_path', type=str, help='Path to output directory where '
                                                                             'extracted neighborhood FASTA files will'
                                                                             ' be saved.')
    return parser.parse_args(args)


def get_filepairs_csv(fasta_path, blast_path, output_path):

    check_output_path(blast_path)

    # Get BLAST results for each AMR gene neighborhood set
    fasta_folder_paths = get_full_filepaths(fasta_path)

    # For each AMR gene, get the neighborhood set for all genomes for a given AMR gene
    header = False
    for AMR_gene_subdir in fasta_folder_paths:

        # Remove non-FASTA files
        remove_files(AMR_gene_subdir, '.fasta')

        # Make output dir for AMR gene within blast output
        AMR_gene_name = get_filename(AMR_gene_subdir)

        # Get all genome neighborhood fasta paths
        root_path = fasta_path + '/' + AMR_gene_name
        full_genome_paths = get_full_filepaths(root_path)

        fields = ['blast_subdir', 'genome_1', 'genome_2']
        row_data = []
        for neighborhood_fasta_1, neighborhood_fasta_2 in itertools.combinations(full_genome_paths, 2):
            row = [blast_path + '/' + AMR_gene_name, neighborhood_fasta_1, neighborhood_fasta_2]
            row_data.append(row)

        with open(output_path + '/genome_pairs.csv', 'a') as csv_file:
            csv_writer = csv.writer(csv_file)
            if not header:
                csv_writer.writerow(fields)
                header = True
            csv_writer.writerows(row_data)


def main(args=None):
    args = parse_args(args)
    get_filepairs_csv(args.FASTA_PATH, args.BLAST_PATH, args.OUTPUT_PATH)


if __name__ == '__main__':
    sys.exit(main())