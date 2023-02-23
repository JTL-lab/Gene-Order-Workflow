#!/usr/bin/env python

"""
Script to generate CSV file containing each unique filepair from among provided whole genome assembly files in 
order to blast them using nf-core DIAMOND within Nextflow workflow.
"""

import argparse
import itertools
import csv
import sys
import os

from utils import get_full_filepaths, remove_files, get_filename, check_output_path


def parse_args(args=None):
    Description = "Create CSV file listing each unique combination of genomes to BLAST in (genome1, genome2) format."
    Epilog = "Example usage: python make_genome_filepairs.py <ASSEMBLY_PATH> <OUTPUT_PATH>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('ASSEMBLY_PATH', metavar='asm_path', type=str,
                        help='Path to directory containing whole genome assembly files to BLAST (e.g. .faa).')
    parser.add_argument('OUTPUT_PATH', metavar='output_path', type=str, help='Path to output directory where '
                                                                             'CSV file will be outputted.')
    return parser.parse_args(args)


def get_filepairs_csv(assembly_path, output_path):

    check_output_path(output_path)

    full_genome_paths = get_full_filepaths(assembly_path)

    fields = ['blast_subdir', 'genome_1', 'genome_2']
    row_data = []
    header = False
    
    # Get every unique combination of genomes to BLAST against each other
    for genome_1, genome_2 in itertools.combinations(full_genome_paths, 2):
        row = [assembly_path, genome_1, genome_2]
        row_data.append(row)

    with open(output_path + '/genome_pairs.csv', 'w') as csv_file:
        csv_writer = csv.writer(csv_file)
        if not header:
            csv_writer.writerow(fields)
            header = True
        csv_writer.writerows(row_data)


def main(args=None):
    args = parse_args(args)
    get_filepairs_csv(args.ASSEMBLY_PATH, args.OUTPUT_PATH)


if __name__ == '__main__':
    sys.exit(main())
