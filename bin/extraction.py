#!/usr/bin/env python

"""
Given RGI files and their corresponding GBK files for complete bacterial genomes,
this script is used for identification of all unique gene neighborhoods present for easy cross-genome comparison.
"""

import os
import glob
import pandas as pd
from Bio import SeqIO
import json
import sys
import argparse
import itertools

from filtering import filter_neighborhoods, write_filtered_genomes_textfile
from utils import get_filename, check_output_path, strip_brackets, generate_alphanumeric_string
from json_utils import make_neighborhood_JSON_data, write_neighborhood_JSON, make_AMR_gene_HTML, \
                       write_clustermap_JSON_HTML

def parse_args(args=None):
    Description = "Extract AMR gene neighborhoods according to fixed window size of N genes upstream and downstream."
    Epilog = "Example usage: python extraction.py <RGI_PATH> <GBK_PATH> <OUTPUT_PATH> -n <NEIGHBORHOOD_SIZE> -p " \
             "<PERCENT> "

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('RGI_PATH', metavar='rgi', type=str, help='Path to directory containing RGI files. Must have '
                                                                  'names corresponding to GBK files.')
    parser.add_argument('GBK_PATH', metavar='gbk', type=str, help='Path to directory containing GBK files. \
                                                                   Must have names corresponding to RGI files.')
    parser.add_argument('OUTPUT_PATH', metavar='output_path', type=str, help='Path to output directory where '
                                                                             'extracted neighborhood FASTA files will'
                                                                             ' be saved.')
    parser.add_argument('-n', metavar='n', type=int, default=10, help='Neighborhood window size, i.e. number of genes '
                                                                      'to consider upstream and downstream of focal '
                                                                      'gene.')
    parser.add_argument('-p', metavar='p', type=float, default=0.75, help='Cutoff percentage of genomes that a given '
                                                                          'AMR gene should be present in for its '
                                                                          'neighborhood to be considered.')
    return parser.parse_args(args)


def load_filepaths(rgi_path_arg, gbk_path_arg):
    """
    Loads all RGI and GBK filepaths from user provided directory paths and returns them in respective lists.
    Assumes that RGI file names have the following naming convention:
    xxxxxxxxxxxxxxxxxxx_rgi.txt
    Assumes that GBK file names have the following naming convention:
    xxxxxxxxxxxxxxxxxxx_genomic.fna.gbk
    """
    # Get paths for RGI files
    try:
        rgi_filepaths = glob.glob(os.path.join(rgi_path_arg, "*.txt"))
    except FileNotFoundError:
        print("Error: there are no RGI files found in the specified directory. Please double-check the provided path.")
        sys.exit(1)

    # Get paths for GBK files
    try:
        gbk_filepaths = glob.glob(os.path.join(gbk_path_arg, "*.gbk"))
    except FileNotFoundError:
        print("Error: there are no GBK files found in the specified directory. Please double-check the provided path.")
        sys.exit(1)

    # Verify there is an equivalent number of RGI and GBK files
    assert len(rgi_filepaths) == len(gbk_filepaths), "Error: mismatch occurred between number of RGI and GBK files."

    # Verify that for each RGI file, there is a GBK file with the same filename
    rgi_without_suffix = [os.path.basename(file).strip('_rgi.txt') for file in rgi_filepaths]
    rgi_list = []
    for file in rgi_without_suffix:
        tokens = file.split('.')
        rgi_list.append(tokens[0])
    rgi_list = [file.split('.')[0] for file in rgi_without_suffix]
    rgi_file_names = set(rgi_list)
    gbk_file_names = set(os.path.basename(file).strip('.gbk') for file in gbk_filepaths)

    # assert rgi_file_names == gbk_file_names, "Error: mismatch occurred between RGI and GBK file names."

    return rgi_filepaths, gbk_filepaths


def load_GBK_file(GBK_filepath):
    """
    Loads a GBK file and gets sequence records, features, and annotations for the file
    """
    filename, extension = os.path.splitext(GBK_filepath)
    assert extension == '.gbk' or extension == '.gb', "Error: filepath provided does not lead to a Genbank file."

    # Extract and store all features and annotations per record
    records = [record for record in SeqIO.parse(GBK_filepath, "genbank")]
    features = [record.features for record in records]
    annotations = [record.annotations for record in records]

    return records, features, annotations


def make_GBK_dataframe(GBK_file_path):
    """
    Input: The output file from parse_genbank_proteins
    Returns a dataframe containing GBK data
    """
    gbk_df = pd.DataFrame()

    gene_start = []
    gene_end = []
    gene_strand = []
    gene_name = []
    loc_tag = []
    function = []
    protein_seq = []
    contig_name = []

    for index, record in enumerate(SeqIO.parse(GBK_file_path, "genbank")):
        for feature in record.features:
            if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers:
                gene_start.append(feature.location.start)
                gene_end.append(feature.location.end)
                gene_strand.append(feature.location.strand)
                loc_tag.append(feature.qualifiers['locus_tag'])
                function.append(feature.qualifiers['product'])
                protein_seq.append(str(feature.qualifiers['translation']))
                contig_name.append(record.id)

                # Preserve gene name if present, otherwise mark as UID (unidentified)
                if 'gene' in feature.qualifiers:
                    gene_name.append(feature.qualifiers['gene'])
                else:
                    gene_name.append("UID-" + str(feature.qualifiers['locus_tag']).strip('[').strip(']'))

    gbk_df['Gene_Start'] = gene_start
    gbk_df['Gene_End'] = gene_end
    gbk_df['Gene_Strand'] = gene_strand
    gbk_df['Locus_Tag'] = loc_tag
    gbk_df['Gene_Name'] = gene_name
    gbk_df['Product'] = function
    gbk_df['Protein_Sequence'] = protein_seq
    gbk_df['Contig_Name'] = contig_name

    unique_contig_names = set(contig_name)

    return gbk_df, unique_contig_names


def make_extraction_dataframe(input_file_path):
    """
    Input: Path to tab-delimited file listing
    Returns a dataframe containing RGI data from RGI
    """
    extraction_df = pd.read_csv(input_file_path, sep='\t', header=0)

    return extraction_df


def adjust_RGI_df_orientation(rgi_df):
    """
    Replaces RGI orientation symbols (-, +) with GBK convention (-1, +1) for easier cross-comparison
    """
    rgi_df['Orientation'] = rgi_df['Orientation'].map({'-': '-1',
                                                       '+': '+1'},
                                                      na_action=None)
    return rgi_df


def shorten_gene_identifier(name):
    """
    Simplifies gene names. Mostly suitable for RGI ARO name cleanup
    (e.g. 'mecC-type BlaZ' -> 'mecC_BlaZ',
          'vanS gene in vanN cluster' -> 'vanS',
          'Bifidobacterium adolescentis rpoB mutants conferring resistance to rifampicin' -> 'rpoB',
          '' -> ''
    )
    """
    conjunctions = ['to', 'and', 'of', 'in']
    tokenized_name = name.split(" ")
    if len(tokenized_name) > 1:
        final_name = tokenized_name[0][0] + tokenized_name[1][0] + "_" + tokenized_name[2]
    else:
        final_name = tokenized_name[0]

    return final_name


def clean_gene_identifier(name):
    """
    Preprocesses gene names for additional shortening/simplification by removing extraneous characters.
    """
    clean_name = name.replace("'", "").replace("-", "_").replace("(", "").replace(")", "").replace("/", "-") \
        .replace(" ", "_").replace(".", "-").replace('"', '').strip('-').strip('"')
    return clean_name


def manipulate_GBK_contigs(df, genome_name):
    """
    For GBK dataframe, manipulates contig to compare if RGI and GBK genes belong to the same contig
    """
    locus_tags = []
    RGI_names = []
    gene_name_identifiers = {}

    df.reset_index(drop=True, inplace=True)
    for index in range(len(df)):
        temp_name = df['Best_Hit_ARO'][index]
        try:
            name = shorten_gene_identifier(temp_name)
        except IndexError:
            name = temp_name
        final_name = clean_gene_identifier(name)
        RGI_names.append(final_name)

        # Keep track of original gene name and the modified version (if any)
        gene_name_identifiers[temp_name] = final_name

        locus_tags.append(genome_name + '(' + final_name + ')' + df['Cut_Off'][index][0] +
                          '_' + str(df['Best_Identities'][index]))

    return locus_tags, RGI_names, gene_name_identifiers


def partition_contig_name(contig_str):
    """
    Splits contig name string and retains first part only (i.e. contig number)
    """
    return contig_str.split("_")[0]


def process_GBK_df_for_BLAST(gbk_df):
    """
    Given a dataframe for a GBK file, removes bracket formatting for easier processing and readability later
    """
    gbk_df['Locus_Tag'] = gbk_df['Locus_Tag'].apply(lambda i: str(i).replace("[", "").replace("]", "").replace("'", ""))
    gbk_df['Gene_Name'] = gbk_df['Gene_Name'].apply(lambda i: str(i).replace("[", "").replace("]", "") \
                                                    .replace("'", "").replace('"', ''))
    gbk_df['ProteinSequence'] = gbk_df['ProteinSequence'].str.strip('[]')

    return gbk_df


def make_RGI_df_contig_col(rgi_df):
    """
    Applies transformation to contig name column to retain only first part of default contig col value from GBK.
    """
    new_contig_col = []
    for value in rgi_df['Contig']:
        head, sep, tail = value.partition("_")
        new_contig_col.append(head)
    rgi_df['Contig'] = new_contig_col

    return rgi_df


def group_by_contig(gbk_df):
    """
    Creates a dict where keys are contig identifiers, and values correspond to all GBK data present for that contig
    """
    group = gbk_df.groupby(gbk_df['Contig_Name'])
    datasets = {}

    for groups, data in group:
        datasets[groups] = data

    return datasets


def get_unique_AMR_genes(RGI_hit_dicts):
    """
    Helper function to return a list of all unique AMR genes present across all inputted genomes.
    """
    unique_AMR_genes = set()
    for RGI_hit_dict in RGI_hit_dicts:
        # Get a list of all AMR genes present in the genome
        AMR_genes = list(RGI_hit_dict.keys())
        # Add to set
        unique_AMR_genes.add(tuple(AMR_genes))

    return unique_AMR_genes


def rename_duplicate_genes(rgi_dataframes):
    """
    Given an RGI dataframe, renames duplicate genes by numbering them according to their contig and start index
    (e.g. multiple instances of vanA become vanA_1-1014, van_2-2500, etc.).
    """
    # Key: genome ID, Val: dict of multiple unique occurrences of the same gene
    multi_gene_genome_dict = {}
    for genome_id, rgi_df in rgi_dataframes.items():

        # Key: gene name, Val = list of contig and start indices for that occurrence
        multi_gene_instances = {}
        unique_genes = []
        for val in range(len(rgi_df)):
            if rgi_df['Best_Hit_ARO'][val] not in unique_genes:
                unique_genes.append(rgi_df['Best_Hit_ARO'][val])

    return


def find_union_AMR_genes(rgi_dataframes):
    """
    TBD with refactoring: used to obtain unique ARO hits from all RGI dataframes and find the union of best hits.
    """
    unique_best_hit_ARO = {}
    union_of_best_hits = []

    # Iterate over every RGI dataframe: need to replace genes in all of them
    for genome_id, rgi_df in rgi_dataframes.items():
        unique_best_hits = []
        for val in range(len(rgi_df)):
            if rgi_df['Best_Hit_ARO'][val] not in unique_best_hits:
                unique_best_hits.append(rgi_df['Best_Hit_ARO'][val])

        unique_best_hit_ARO[genome_id] = list(set(unique_best_hits))

    for genome_id, ARO_term in unique_best_hit_ARO.items():
        for val in ARO_term:
            if val not in union_of_best_hits:
                union_of_best_hits.append(val)

    return unique_best_hit_ARO, union_of_best_hits


def check_duplicate_genes(RGI_dataframes, ARO_union):
    """
    Given the RGI dataframes for the genomes being analyzed, labels multi-gene instances numerically to differentiate
    between them downstream.
    """
    all_unique_ARO = []
    for AMR_gene in ARO_union:
        for genome, rgi_df in RGI_dataframes.items():

            # Locate all instances of the gene within the genome
            AMR_gene_row = rgi_df.loc[rgi_df['Best_Hit_ARO'] == AMR_gene]

            # Case 1: If multiple instances are present, make a new dictionary for each new gene and modify the df
            if len(AMR_gene_row) > 1:
                instance_num = 1
                for gene_instance in range(len(AMR_gene_row)):
                    # Relabel each instance and save in the df
                    gene_name_tokens = rgi_df['Best_Hit_ARO'][gene_instance].split('_')
                    if instance_num > 1:
                        del gene_name_tokens[-1]
                    restored_gene_name = '_'.join(gene_name_tokens)
                    gene_name = restored_gene_name + '_' + str(instance_num)
                    rgi_df.loc[:, ('Best_Hit_ARO', gene_instance)] = gene_name
                    instance_num += 1

                    # Append with new name
                    all_unique_ARO.append(gene_name)

            elif len(AMR_gene_row) == 1:
                all_unique_ARO.append(AMR_gene)

    return RGI_dataframes, set(all_unique_ARO)


def make_AMR_dict(RGI_dataframes, AMR_gene_index):
    """
    Given the RGI dataframes for the genomes being analyzed, creates a dictionary of AMR genes
    """
    AMR_dict = {}
    for genome, rgi_df in RGI_dataframes.items():
        AMR_gene_row = rgi_df.loc[rgi_df['Best_Hit_ARO'] == AMR_gene_index]
        if len(AMR_gene_row) > 0:
            AMR_dict[genome] = AMR_gene_row

    return AMR_dict


def make_AMR_gene_neighborhood_df(GBK_df_dict, genome_id, gene_start, gene_name, neighborhood_size):
    """
    Finds gene neighborhood of size 2N (user-defined by neighborhood_size, i.e., N genes downstream and upstream)
    for a given AMR gene for cross genome comparison.
    """
    # Get the GBK data for the given genome
    try:
        GBK_df = GBK_df_dict[genome_id]
        GBK_df.reset_index(drop=True, inplace=True)
    except KeyError:
        return TypeError

    # Keep track of contig ends: track relative start and stop coordinates of neighborhood
    neighborhood_indices = []

    # Get the focal (central) AMR gene
    try:
        # Subtract one from gene start index to account for automatic padding
        AMR_gene_df_row = GBK_df.loc[((GBK_df['Gene_Start'] == gene_start - 1))]

        # Get only genes on the same contig as the focal gene for consideration as neighbors
        contig_id = AMR_gene_df_row.Contig_Name.tolist()[0]
        contig_df = GBK_df.loc[(GBK_df['Contig_Name'] == contig_id)].sort_values(by='Gene_Start')
        contig_df.reset_index(drop=True, inplace=True)

        AMR_gene_index = contig_df.index[(contig_df['Gene_Start'] == gene_start - 1)].tolist()
        gene_index = AMR_gene_index[0]

        # Get downstream neighbors
        downstream = [gene_index - index for index in range(1, neighborhood_size + 1)]

        # If contig end present downstream, some indices will be negative: remove these to prevent index errors
        downstream_indices = [index for index in downstream if index >= 0]
        downstream_neighbors = pd.DataFrame(columns=['Gene_Start', 'Gene_End', 'Gene_Strand', 'Locus_Tag', 'Gene_Name',
                                                     'Product', 'Protein_Sequence', 'Contig_Name'])
        for i in range(len(downstream_indices) - 1, -1, -1):
            try:
                neighbor = contig_df.iloc[[downstream_indices[i]]]
                downstream_neighbors = pd.concat([downstream_neighbors, neighbor])
            except IndexError:
                print("Contig end found at position -{} downstream.".format(i + 1))
                neighborhood_indices.append(i + 1)

        # If there was no contig end, append default N size to neighborhood_indices
        if len(neighborhood_indices) == 0:
            neighborhood_indices.append(-neighborhood_size)

        # Get upstream neighbors
        upstream_indices = [gene_index + index for index in range(1, neighborhood_size + 1)]
        upstream_neighbors = pd.DataFrame(columns=['Gene_Start', 'Gene_End', 'Gene_Strand', 'Locus_Tag',
                                                   'Gene_Name', 'Product', 'Protein_Sequence', 'Contig_Name'])
        for i in range(len(upstream_indices)):
            contig_end_found = False
            try:
                neighbor = contig_df.iloc[[upstream_indices[i]]]
                upstream_neighbors = pd.concat([upstream_neighbors, neighbor])
            except IndexError:
                if not contig_end_found:
                    print("Contig end found at position {} upstream.".format(i + 1))
                    contig_end_found = True
                if len(neighborhood_indices) < 2:
                    neighborhood_indices.append(i + 1)

        if len(neighborhood_indices) < 2:
            neighborhood_indices.append(neighborhood_size)

        neighborhood_df = pd.concat([upstream_neighbors, AMR_gene_df_row, downstream_neighbors])

        return neighborhood_df, neighborhood_indices

    except IndexError:
        print("Gene {gene} not found.".format(gene=gene_name))


def get_all_AMR_gene_neighborhoods(AMR_instance_dict, GBK_df_dict, unique_AMR_genes, neighborhood_size):
    """
    Given a dictionary of AMR genes determined using RGI outputs and a list of unique AMR genes, creates one dictionary
    containing all AMR gene neighborhoods for a fixed size of 2N (user-defined by neighborhood_size, i.e., N genes
    downstream and upstream) for the set of genomes being analyzed.
    """
    # Will be used to store dataframes of neighborhoods
    gene_neighborhoods = {}

    # Keeps track of contig ends, if applicable
    contig_ends = {}

    for AMR_gene, AMR_dict in AMR_instance_dict.items():

        # Keep track of size 2N neighborhood
        neighbors = {}

        # Get the name without any numbering for multiple instances, since start index will allow us to differentiate
        gene_name = AMR_gene.split('_')[0]

        # Keep track of start and stop indices for each neighborhood
        contig_end_flags = {}
        errors = 0

        for genome, data in AMR_dict.items():

            # Get start index of the focal AMR gene from the RGI dataframe
            start_vals = AMR_dict[genome]['Start']
            start_vals_list = list(start_vals)
            start_index = start_vals_list[0]

            # Make gene neighborhood dataframe for each genome for the focal gene, AMR_gene
            try:
                neighborhood_df, indices = make_AMR_gene_neighborhood_df(GBK_df_dict, genome, start_index, gene_name,
                                                                         neighborhood_size)
                neighborhood_df.reset_index(drop=True, inplace=True)

                if len(neighborhood_df) > 1:
                    neighbors[genome] = neighborhood_df
                    contig_end_flags[genome] = indices
            except TypeError:
                errors += 1
                print("AMR gene {} was not present in this genome!".format(AMR_gene))

        gene_neighborhoods[AMR_gene] = neighbors
        contig_ends[AMR_gene] = contig_end_flags

    return gene_neighborhoods, contig_ends


def get_neighborhood_gene_data(neighborhood_df):
    """
    Given a neighborhood dataframe, obtains the locus tags, protein sequences, and gene names as separate dicts
    """
    locus_list = neighborhood_df['Locus_Tag'].tolist()
    protein_list = neighborhood_df['Protein_Sequence'].tolist()
    gene_name_list = neighborhood_df['Gene_Name'].tolist()

    locus_to_protein_dict = {}
    for gene in range(len(locus_list)):
        locus_to_protein_dict[locus_list[gene].strip()] = protein_list[gene]

    return locus_list, protein_list, gene_name_list


def get_neighborhood_data(neighborhoods_dict, num_neighbors):
    """
    Extracts locus, protein sequence, and gene name data for a dictionary of neighborhood data created using
    get_all_AMR_gene_neighborhoods.
    """
    locus_dict = {}
    protein_dict = {}
    gene_name_dict = {}

    for AMR_gene, genome_neighborhoods in neighborhoods_dict.items():
        locus_tags = {}
        protein_seqs = {}
        gene_names = {}

        for genome_id, neighborhood_df in genome_neighborhoods.items():
            locus_tags[genome_id], protein_seqs[genome_id], gene_names[genome_id] = \
                get_neighborhood_gene_data(genome_neighborhoods[genome_id])

        locus_dict[AMR_gene] = locus_tags
        protein_dict[AMR_gene] = protein_seqs
        gene_name_dict[AMR_gene] = gene_names

    return locus_dict, protein_dict, gene_name_dict


def write_AMR_neighborhood_to_FNA(AMR_gene_neighborhoods_dict, AMR_gene, locuses_dict, protein_seqs_dict, out_path):
    """
    Given a dictionary containing AMR gene neighborhoods for a set of genomes being analyzed and a specified output
    directory path, creates a distinct .fna file for each AMR gene neighborhood to use for BLAST All-vs-All comparison.
    """
    for genome_id, genome_neighborhood_df in AMR_gene_neighborhoods_dict.items():
        # Create a new FASTA file to write neighborhood sequence and identifier data to
        with open(out_path + '/' + genome_id + '.fasta', 'w') as output_fasta:
            for locus, protein_seq in zip(locuses_dict[AMR_gene][genome_id], protein_seqs_dict[AMR_gene][genome_id]):
                locus_output = locus.strip("'")
                protein_output = protein_seq.strip("'")
                output_fasta.write('>{}\n'.format(locus_output))
                output_fasta.write('{}\n'.format(protein_output))
        output_fasta.close()


def delete_low_occurring_genes(AMR_gene_dict, num_genomes, cutoff_percentage=0.25):
    """
    Removes keys of AMR genes not present in cutoff percentage of genomes.
    - AMR_gene_dict should be the dictionary of one type of AMR gene, with entries representing
    its occurrences within genomes.
    - Cutoff percentage is a float between 0 and 1 (e.g., 0.25 means only genes present in min 25% of genomes are kept).
    """
    # Define threshold amount of genomes an AMR gene must be present in for inclusion
    minimum_num_genomes = round(num_genomes * cutoff_percentage)

    # Remove all AMR genes that occur in less genomes than the threshold amount
    genes_to_remove = []
    for AMR_gene, neighborhoods_dict in AMR_gene_dict.items():
        if len(neighborhoods_dict) < minimum_num_genomes:
            genes_to_remove.append(AMR_gene)

    for AMR_gene in genes_to_remove:
        del AMR_gene_dict[AMR_gene]

    return AMR_gene_dict


def get_AMR_gene_statistics(rgi_dataframes):
    """
    Given RGI dataframes, counts the number of Perfect and Strict hits in each genome.
    Needed for extraction summary file.
    """
    amr_gene_statistics = {}
    for genome, rgi_df in rgi_dataframes.items():
        try:
            perfect_count = len(rgi_df[rgi_df["Cut_Off"] == "Perfect"])
            strict_count = len(rgi_df[rgi_df["Cut_Off"] == "Strict"])
            loose_count = len(rgi_df[rgi_df["Cut_Off"] == "Loose"])
            amr_gene_statistics[genome] = [perfect_count, strict_count, loose_count]
        except IndexError:
            pass

    return amr_gene_statistics


def write_summary_file(output_dir, gbk_dataframes, neighborhood_size, ARO_union,
                       amr_gene_statistics, gene_name_identifiers_dict):
    """
    Writes summary data of neighborhood extraction for the set of genomes being analyzed to
    output/Neighborhoods_Extraction_Summary.txt.
    """
    with open(output_dir + '/' + 'Neighborhoods_Extraction_Summary.txt', 'w') as outfile:
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        outfile.write('NEIGHBORHOOD EXTRACTION SUMMARY' + '\n')
        outfile.write(
            '----------------------------------------------------------------------------------------' + '\n\n')
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        outfile.write('General Details' + '\n')
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        outfile.write('Number of genomes analyzed: {}'.format(len(gbk_dataframes.keys())) + '\n')
        outfile.write('Neighborhood size extracted: {}'.format(neighborhood_size) + '\n')
        outfile.write('Number of genes whose neighborhoods were extracted: {}'.format(len(ARO_union)) + '\n\n')
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        outfile.write('Original Gene Names and Simplified Names' + '\n')
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        for gene_name, simplified_gene_name in sorted(gene_name_identifiers_dict.items()):
            outfile.write(gene_name + '\t--->\t' + simplified_gene_name + '\n')
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        outfile.write('AMR Gene Details Per Genome' + '\n')
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        outfile.write('GENOME_ID' + '\t' + 'PERFECT_HITS' + '\t' + 'STRICT_HITS' + '\t' + 'LOOSE_HITS' + '\n')
        for genome, counts in amr_gene_statistics.items():
            outfile.write('{genome}\t{perfect}\t{strict}\t{loose}\n'.format(genome=get_filename(genome),
                                                                            perfect=counts[0],
                                                                            strict=counts[1],
                                                                            loose=counts[2]))


def extract_neighborhoods(rgi_path, gbk_path, output_path, num_neighbors, cutoff_percent):
    """
    Driver script for extracting all AMR gene neighborhoods from specified GBK and RGI files to output FASTA format
    protein sequences for each neighborhood.
    """
    # 0) Get user path args from pipeline for path to RGI files and GBK files
    try:
        rgi_filepaths, gbk_filepaths = load_filepaths(rgi_path, gbk_path)
    except IndexError:
        print("Please provide paths to the RGI files and GBK files respectively when running this script.")
        sys.exit(1)

    print("RGI filepaths: {}".format(rgi_filepaths))
    print("GBK filepaths: {}".format(gbk_filepaths))

    # Make output directory if non-existent
    check_output_path(output_path)

    print("Processing GBK and RGI files...")

    # 1) Get Genbank records, features and annotations for each GBK file
    # Contains a dataframe for each GBK file, where keys are filenames
    gbk_dataframes = {}

    # Contains contig names associated with each GBK file required for later processing, where keys are filenames
    gbk_contigs = {}

    for gbk_filepath in gbk_filepaths:

        # Get GBK filename for GBK dictionary keys
        gbk_filename = get_filename(gbk_filepath)

        # Get GBK dataframe and contig names
        gbk_df, contig_names = make_GBK_dataframe(gbk_filepath)

        # Preprocess all bracketed columns to remove brackets
        for col in gbk_df:
            if type(gbk_df[col][0]) is list or isinstance(gbk_df[col][0], str):
                gbk_df[col] = gbk_df[col].apply(lambda x: strip_brackets(x))

        # Retain GBK dataframe and GBK list of contig names with filename as key
        gbk_dataframes[gbk_filename] = gbk_df
        gbk_contigs[gbk_filename] = contig_names

    # 2) Get RGI dataframes and do all required preprocessing for later manipulation and comparison
    # Contains a dataframe for each RGI file, where keys are filenames
    rgi_dataframes = {}

    for rgi_filepath in rgi_filepaths:
        rgi_filename = get_filename(rgi_filepath, rgi=True)
        rgi_df = make_extraction_dataframe(rgi_filepath)
        rgi_dataframes[rgi_filename] = rgi_df

    for rgi_filename, rgi_df in rgi_dataframes.items():

        # Replace each RGI dataframe with one with the proper DNA orientation representation to match GBKs
        rgi_dataframes[rgi_filename] = adjust_RGI_df_orientation(rgi_df)

        # Add locus tag and contig details to the RGI dataframes
        rgi_dataframes[rgi_filename] = make_RGI_df_contig_col(rgi_df)
        rgi_dataframes[rgi_filename]['Locus_Tag'], rgi_dataframes[rgi_filename][
            'Best_Hit_ARO'], gene_name_identifiers_dict = manipulate_GBK_contigs(rgi_df, rgi_filename)

        # Preprocess all bracketed columns to remove brackets
        for col in rgi_df:
            if type(rgi_df[col][0]) is list or isinstance(rgi_df[col][0], str):
                try:
                    rgi_df[col] = rgi_df[col].apply(lambda x: strip_brackets(x))
                except AttributeError:
                    pass

    # 4) Get unique ARO terms and union of all unique ARO terms
    ARO_best_terms, ARO_union = find_union_AMR_genes(rgi_dataframes)

    # 5) Get AMR dict for all genomes
    # rgi_modified_dataframes, unique_ARO = check_duplicate_genes(rgi_dataframes, ARO_union)
    # genomes_AMR_dict = make_AMR_dict(rgi_modified_dataframes, unique_ARO)
    genomes_AMR_dict = {}
    for AMR_gene in ARO_union:
        # Make entry for gene occurrences
        genomes_AMR_dict[AMR_gene] = make_AMR_dict(rgi_dataframes, AMR_gene)

    amr_statistics = get_AMR_gene_statistics(rgi_dataframes)

    genomes_AMR_dict_filtered = delete_low_occurring_genes(genomes_AMR_dict, len(gbk_filepaths), cutoff_percent)

    print("Extracting gene neighborhood data for neighborhoods of size {}...".format(num_neighbors))

    # 6) Extract AMR neighborhoods and store them in neighborhood dataframes
    neighborhoods, neighborhood_indices = get_all_AMR_gene_neighborhoods(genomes_AMR_dict_filtered, gbk_dataframes,
                                                                         ARO_union, num_neighbors)

    # 7) Get the locus, protein, and gene name details for each neighborhood respectively for FNA file creation
    locuses, protein_seqs, gene_names = get_neighborhood_data(neighborhoods, num_neighbors)

    # 8) Save neighborhoods to FNA files needed for All-vs-all BLAST results later
    print("Generating neighborhood FNA files...")
    for AMR_gene, genome_neighborhoods in neighborhoods.items():

        # Make output subdirectory for the gene
        gene_name = AMR_gene.replace('(', '_').replace(')', '_').replace("'", "")
        out_path = output_path + '/fasta/' + gene_name
        check_output_path(out_path)

        for genome_id in genome_neighborhoods.keys():
            # Make file with genome ID as filename
            write_AMR_neighborhood_to_FNA(genome_neighborhoods, AMR_gene, locuses, protein_seqs, out_path)

    # 9) Make neighborhoods summary textfile in output dir
    print("Making extraction summary file...")
    write_summary_file(output_path, gbk_dataframes, num_neighbors, ARO_union, amr_statistics, gene_name_identifiers_dict)

    # 10) Save gene neighborhoods and indices in textfile: needed for JSON representations downstream
    for AMR_gene, neighborhood_data in neighborhoods.items():
        # Filter neighborhoods
        filtered_neighborhoods = filter_neighborhoods(neighborhood_data, num_neighbors)
        surrogates = list(filtered_neighborhoods.keys())
        filtered_neighborhoods_dict = {key: value for (key, value) in neighborhood_data.items() if key in surrogates}
        write_filtered_genomes_textfile(filtered_neighborhoods, AMR_gene, output_path)

        # Get data needed to write JSON files
        neighborhood_JSON_dict = make_neighborhood_JSON_data(neighborhood_data, AMR_gene)
        filtered_neighborhood_JSON_dict = make_neighborhood_JSON_data(filtered_neighborhoods_dict, AMR_gene)

        # Update UID genes' data if present
        # complete_neighborhood_JSON_dict = update_UID_JSON_data(neighborhood_JSON_dict, AMR_gene, output_path)

        # Create JSON file
        # write_neighborhood_JSON(complete_neighborhood_JSON_dict, AMR_gene, output_path)
        write_neighborhood_JSON(neighborhood_JSON_dict, AMR_gene, output_path)
        write_neighborhood_JSON(filtered_neighborhood_JSON_dict, AMR_gene, output_path, True)

    sample_data_path = '../../../sample_data'
    make_AMR_gene_HTML(neighborhoods.keys(), sample_data_path, output_path)

    with open(output_path + '/' + 'neighborhood_indices.txt', 'w') as outfile:
        outfile.write(str(neighborhood_indices))

    print("Neighborhood extraction complete.")
    check_output_path(output_path + '/blast')


def main(args=None):
    args = parse_args(args)
    extract_neighborhoods(args.RGI_PATH, args.GBK_PATH, args.OUTPUT_PATH, args.n, args.p)


if __name__ == '__main__':
    sys.exit(main())
