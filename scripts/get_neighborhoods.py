#!/usr/bin/env python

"""
Given RGI files and their corresponding GBK files for complete bacterial genomes,
this script is used for identification of all unique gene neighborhoods present for easy cross-genome comparison.

Acknowledgement: This script is largely adapted from the work of C. N. Rudrappa, which can be found at:
https://github.com/chandana277/AMR_Analysis_Pipeline/blob/main/Neighbor_Generator.py
"""

import os
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import json
import re
import sys
import warnings
import argparse
import time

def load_filepaths(rgi_path_arg, gbk_path_arg):
    """
    Loads all RGI and GBK filepaths from user provided directory paths and returns them in respective lists.

    Assumes that RGI file names have the following naming convention:
    xxxxxxxxxxxxxxxxxxx_genomic.fna_rgi.txt

    Assumes that GBK file names have the following naming convention:
    xxxxxxxxxxxxxxxxxxx_genomic.fna.gbk
    """

    # Get paths for RGI files
    try:
        rgi_filepaths = glob.glob(os.path.join(rgi_path_arg, "*.txt"))
    except FileNotFoundError:
        print("Error: there are no RGI files found in the specified directory. Please double-check the provided path.")

    # Get paths for GBK files
    try:
        gbk_filepaths = glob.glob(os.path.join(gbk_path_arg, "*.gbk"))
    except FileNotFoundError:
        print("Error: there are no GBK files found in the specified directory. Please double-check the provided path.")

    # Verify there is an equivalent number of RGI and GBK files
    assert len(rgi_filepaths) == len(gbk_filepaths), "Error: mismatch occurred between number of RGI and GBK files."

    # Verify that for each RGI file, there is a GBK file with the same filename
    rgi_file_names = set(os.path.basename(file).strip('.fna_rgi.txt') for file in rgi_filepaths)
    gbk_file_names = set(os.path.basename(file).strip('.gbk') for file in gbk_filepaths)

    assert rgi_file_names == gbk_file_names, "Error: mismatch occurred between RGI and GBK file names."

    return rgi_filepaths, gbk_filepaths

def get_filename(filepath):
    """
    Gets filename from filepath (i.e. removes directories and file suffix)
    """
    filename = os.path.basename(filepath).split(".")[0]
    return filename

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

def make_GBK_dict(GBK_filepath):
    """
    Makes a dictionary of gene names with their start index as the key
    """
    GBK_dict = {}
    records, features, annotations = load_GBK_file(GBK_filepath)

    for record in records:
        for feature in record.features:
            if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers:

                # Get start index of gene
                location = str(feature.location).split(':')
                start = location[0].strip('[')

                # Get gene name
                gene_name = feature.qualifiers['gene']

                GBK_dict[start] = gene_name

    return GBK_dict

def parse_genbank_proteins(features_list):
    """
    Obtains proteins from GBK file features_list
    """
    protein_features = []
    for feature in features_list:
        # Keep only valid proteins!
        if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers:
            protein_features.append(feature)

    return protein_features

def get_protein_sequences(records, protein_features):
    """
    Makes a list of protein sequences from a list of protein features
    """
    seqs = []
    for record in records:
        for feature in protein_features:
            seqs.append(SeqRecord(Seq(record.qualifiers['translation'][0]),
                                  id=record.qualifiers['locus_tag'][0],
                                  name=record.name))
    return seqs
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
                    gene_name.append("UID")

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

def make_RGI_dataframe(rgi_filepath):
    """
    Input: path to RGI file
    Returns a dataframe containing RGI data from RGI
    """
    rgi_df = pd.read_csv(rgi_filepath, sep='\t', header=0)

    return rgi_df

def strip_brackets(var):
    """
    Preprocesses list or string with brackets to be bracketless string for easier manipulation of dataframes
    """
    if isinstance(var, list):
        final = str(var)
        clean_str = final.strip('[]')
    else:
        clean_str = var.strip('[]')

    return clean_str

def adjust_RGI_df_orientation(rgi_df):
    """
    Replaces RGI orientation symbols (-, +) with GBK convention (-1, +1) for easier cross-comparison
    """
    rgi_df['Orientation'] = rgi_df['Orientation'].map({'-': '-1',
                                                       '+': '+1'},
                                                      na_action=None)
    return rgi_df

def manipulate_GBK_contigs(dataframe, genome_name):
    """
    For GBK dataframe, manipulates contig to compare if RGI and GBK genes belong to the same contig
    """
    locus_tags = []
    RGI_names = []

    dataframe.reset_index(drop=True, inplace=True)
    dataframe.reset_index(drop=True, inplace=True)
    for index in range(len(dataframe)):
        temp_name = dataframe['Best_Hit_ARO'][index].split(" ")

        if len(temp_name) > 1:
            name = temp_name[0][0] + temp_name[1][0] + "_" + temp_name[2]
            RGI_names.append(name)
        else:
            name = temp_name[0]
            RGI_names.append(name)

        locus_tags.append(genome_name + '(' + name + ')' + dataframe['Cut_Off'][index][0] + '_' + str(dataframe['Best_Identities'][index]))

    return locus_tags, RGI_names

def partition_contig_name(contig_str):
    """
    Splits contig name string and retains first part only (i.e. contig number)
    """
    return contig_str.split("_")[0]

def process_GBK_df_for_BLAST(gbk_df):
    """
    Given a dataframe for a GBK file, removes
    """
    gbk_df['Locus_Tag'] = gbk_df['Locus_Tag'].apply(lambda i:str(i).replace("[","").replace("]","").replace("'",""))
    gbk_df['ProteinSequence'] = gbk_df['ProteinSequence'].str.strip('[]')

    return gbk_df

def make_RGI_df_contig_col(rgi_df):
    """
    Applies transformation to contig name column to retain only first part of default contig col value from GBK
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

def find_union_AMR_genes(rgi_dataframes):
    """
    TBD with refactoring: used to obtain unique ARO hits from all RGI dataframes and find the union of best hits
    """
    unique_best_hit_ARO = {}
    union_of_best_hits = []

    for genome_id, rgi_df in rgi_dataframes.items():
        unique_best_hits = []
        for val in range(len(rgi_df)):
            if rgi_df['Best_Hit_ARO'][val] not in unique_best_hits:
                unique_best_hits.append(rgi_df['Best_Hit_ARO'][val])
        unique_best_hit_ARO[genome_id] = unique_best_hits

    for genome_id, ARO_term in unique_best_hit_ARO.items():
        for val in ARO_term:
            if val not in union_of_best_hits:
                union_of_best_hits.append(val)

    return unique_best_hit_ARO, union_of_best_hits
def get_RGI_instances(AMR_dicts):
    """
    Get single instance and multiple instance (i.e. occurring in multiple places in genome data) AMR genes respectively
    """
    single_instance_AMR = {}
    multiple_instance_AMR = {}

    # Check the AMR dictionary for each respective AMR gene found
    for AMR_gene, AMR_dict in AMR_dicts.items():
        for key, val in AMR_dict.items():
            # If AMR gene occurred in multiple places, add to multiple instance dict to keep track of this
            if len(val) > 1:
                multiple_instance_AMR[AMR_gene] = AMR_dict
            else:
                single_instance_AMR[AMR_gene] = AMR_dict

    return single_instance_AMR, multiple_instance_AMR
def get_neighborhood_gene_data(neighborhood_df, num_genes):
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
    GBK_df = GBK_df_dict[genome_id]
    GBK_df.reset_index(drop=True, inplace=True)

    contig_flag = 0
    # Get the focal (central) AMR gene
    try:
        print("GENE START: {}".format(gene_start))
        # Subtract one from gene start index to account for automatic padding
        AMR_gene_df_row = GBK_df.loc[(GBK_df['Gene_Start'] == gene_start - 1)]
        AMR_gene_index = GBK_df.index[(GBK_df['Gene_Start'] == gene_start - 1)].tolist()
        gene_index = AMR_gene_index[0]

        # Get downstream neighbors
        downstream_indices = [gene_index - index for index in range(1, neighborhood_size + 1)]
        downstream_neighbours = GBK_df.iloc[downstream_indices]

        # Get upstream neighbors
        upstream_indices = [gene_index + index for index in range(1, neighborhood_size + 1)]
        upstream_neighbours = GBK_df.iloc[upstream_indices]

        # Create a dataframe containing the full neighborhood data; TBD verify that contig ends handled
        neighborhood_df = pd.concat([downstream_neighbours, AMR_gene_df_row, upstream_neighbours]).sort_values(by="Gene_Start")

        return neighborhood_df

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
        neighbours = {}

        # Keep track of contig ends for a given neighborhood
        contig_end_flags = {}

        for genome, data in AMR_dict.items():

            # Get start index of the focal AMR gene from the RGI dataframe
            start_vals = AMR_dict[genome]['Start']

            start_vals_list = list(start_vals)
            start_index = start_vals_list[0]

            # Make gene neighborhood dataframe for each genome for the focal gene, AMR_gene
            neighborhood_df = make_AMR_gene_neighborhood_df(GBK_df_dict, genome, start_index,
                                                            AMR_gene, neighborhood_size)

            try:
                if len(neighborhood_df) > 0:
                    neighbours[genome] = neighborhood_df
            except TypeError:
                print("AMR gene {} was not present in this genome!".format(AMR_gene))

        gene_neighborhoods[AMR_gene] = neighbours
        contig_ends[AMR_gene] = contig_end_flags

    return gene_neighborhoods, contig_ends

def get_neighborhood_data(instance_neighborhoods_dict):
    """
    Extracts locus, protein sequence, and gene name data for a dictionary of neighborhood data created using
    get_all_AMR_gene_neighborhoods.
    """
    instance_locus_data = {}
    instance_protein_data = {}
    instance_gene_name_data = {}
    for neighborhood, neighborhood_data in instance_neighborhoods_dict.items():
        locus_tags = {}
        protein_seqs = {}
        gene_names = {}

        for gene, neighbor in neighborhood_data.items():
            locus_tags[gene], protein_seqs[gene], gene_names[gene] = get_neighborhood_gene_data(neighborhood_data[gene],
                                                                                                num_neighbours)

        instance_locus_data[neighborhood] = locus_tags
        instance_protein_data[neighborhood] = protein_seqs
        instance_gene_name_data[neighborhood] = gene_names

    return instance_locus_data, instance_protein_data, instance_gene_name_data

def write_AMR_neighborhood_to_FNA(AMR_neighborhoods_dict, AMR_gene, genome_id, locuses_dict,
                                  protein_sequences_dict, output_path):
    """
    Given an dictionary containing AMR gene neighborhoods for a set of genomes being analyzed and a specified output
    directory path, creates a distinct .fna file for each AMR gene neighborhood to use for BLAST All-vs-All comparison.
    """
    neighborhood_dict = AMR_neighborhoods_dict[AMR_gene]
    # Create a new FASTA file to write neighborhood sequence and identifier data to
    with open(output_path + '/' + genome_id + '.fasta', 'w') as output_fasta:
        for neighbor in neighborhood_dict.keys():
            for locus, protein_seq in zip(locuses_dict[AMR_gene][neighbor], protein_sequences_dict[AMR_gene][neighbor]):
                output_fasta.write('>{}\n'.format(locus))
                output_fasta.write('{}\n'.format(protein_seq))
        output_fasta.close()

def get_RGI_hits(rgi_filepath, gbk_dict):
    """
    Returns a dictionary of all RGI hits for a genome, along with a list of their contig start and stop indices
    """
    AMR_dict = {}

    # Iterate over RGI file lines (skipping the header)
    with open(rgi_filepath, 'r') as rgi_input:
        for line in itertools.islice(rgi_input, 1, None):

            AMR_hit = []

            rgi_tokens = re.split(" |\t|;", line)

            # Add hit type, start, and stop indices for the AMR gene
            if 'Perfect' in rgi_tokens:
                AMR_hit.append('Perfect')
            elif 'Strict' in rgi_tokens:
                AMR_hit.append('Strict')
            else:
                AMR_hit.append('Loose')

            # Add start and stop indices for the feature
            start = rgi_tokens[2]
            stop = rgi_tokens[4]

            AMR_hit.append(start)
            AMR_hit.append(stop)

            # Get gene name
            gene_name = gbk_dict[start]

            # Retain in dictionary
            AMR_dict[gene_name] = AMR_hit

    rgi_input.close()

    return AMR_dict

def make_AMR_presence_dict(RGI_hit_dicts, unique_AMR_genes):
    """
    Given a list of RGI hit dicts for each genome created using the get_RGI_hits function AND a set of the unique AMR
    genes present in all genomes, returns a dictionary where AMR genes are keys and values are a list of all genomes
    they are present in.
    """
    AMR_presence_dict = {}
    for AMR_gene in unique_AMR_genes:
        for RGI_dict in RGI_hit_dicts:

            # Try looking up the gene: if present, retain the genome name, otherwise ignore
            try:
                indices = RGI_dict[AMR_gene]
            except IndexError:
                pass

    return AMR_presence_dict

def delete_low_occurring_genes(AMR_gene_dict, num_genomes, cutoff_percentage=0.25):
    """
    Removes keys of AMR genes not present in cutoff percentage of genomes (recommended minimum/default of 25%).
    - AMR_gene_dict should be the dictionary of one type of AMR gene, with entries representing
    its occurrences within genomes.
    - Cutoff percentage is a float between 0 and 1 (e.g., 0.8 means only genes present in min 80% of genomes are kept).
    """
    # Define threshold amount of genomes an AMR gene must be present in for inclusion
    minimum_num_genomes = num_genomes * cutoff_percentage

    # Remove all AMR genes that occur in less genomes than the threshold amount
    for item, occurrences in AMR_gene_dict.items():
        if len(occurrences) < minimum_num_genomes:
            del AMR_gene_dict[item]

    return AMR_gene_dict

def make_gene_neighborhood_JSON(AMR_gene_neighborhood_sets):
    """
    Creates JSON file containing neighborhood data for an AMR gene for the genomes under analysis of the following form
    (where the first N neighbour genes represent the N neighbours downstream of the target gene, and the last N
    neighbour genes represent the N neighbours upstream of the target gene):

    { "gene_name" : "geneA",
      "neighbour_genes": [
        "genome_1" : ["geneB", "geneC", ..., "geneZ"]
        "genome_2" : ["geneC", "geneB", ..., "geneY"]
                              ...
        "genome_X" : ["geneY", "geneB", ..., "geneC"]
      ],
     "neighbour_indices": [
        "genome_1" : [[1000, 1238], [1238, 1334], ..., [2345, 2455]]
        "genome_2" : [[1100, 1248], [1278, 1374], ..., [2445, 2555]]
                              ...
        "genome_X" : [[1100, 1248], [1278, 1374], ..., [2445, 2555]]
     ]
    }

    This function should be run on the complete set of neighborhoods (i.e. for all AMR genes).
    """

    with open('neighborhoods.json', 'a') as outfile:
        json.dump(AMR_gene_neighborhood_sets, outfile)
    outfile.close()
    return

def check_output_path(path):
    """
    Checks for presence of specified output path and creates the directory if it does not exist
    """
    if not os.path.exists(path):
        os.makedirs(path)

def write_summary_file(output_dir, num_genomes, neighborhood_size):
    """
    Writes summary data of neighborhood extraction for the set of genomes being analyzed to
    output/Neighborhoods_Extraction_Summary.txt.
    """
    with open(output_dir + '/' + 'Neighborhoods_Extraction_Summary.txt', 'w') as outfile:
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        outfile.write('NEIGHBORHOODS SUMMARY' + '\n')
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        outfile.write('Number of genomes analyzed: {}'.format(num_genomes) + '\n')
        outfile.write('Neighborhood size extracted: {}'.format(neighborhood_size) + '\n')
    outfile.close()

if __name__ == '__main__':

    start_time = time.time()

    # 0) Get user path args from pipeline for path to RGI files and GBK files
    try:
        rgi_path = sys.argv[1]
        gbk_path = sys.argv[2]
    except IndexError:
        print("Please provide paths to the RGI files and GBK files respectively when running this script.")

    rgi_filepaths, gbk_filepaths = load_filepaths('sample_data/RGI', 'sample_data/GBK')

    try:
        num_neighbours = int(sys.argv[3])
    except IndexError:
        # Default value
        num_neighbours = 10

    try:
        cutoff_percent = float(sys.argv[4])
    except IndexError:
        # Provide default value
        cutoff_percent = 0.25

    # 1) Get Genbank records, features and annotations for each GBK file
    records = []
    features = []
    annotations = []

    for gbk_file in gbk_filepaths:
        rec, feat, annot = load_GBK_file(gbk_file)
        records.append(rec)
        features.append(feat)
        annotations.append(annot)

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
        rgi_filename = get_filename(rgi_filepath)
        rgi_df = make_RGI_dataframe(rgi_filepath)
        rgi_dataframes[rgi_filename] = rgi_df

    for rgi_filename, rgi_df in rgi_dataframes.items():

        # Replace each RGI dataframe with one with the proper DNA orientation representation to match GBKs
        rgi_dataframes[rgi_filename] = adjust_RGI_df_orientation(rgi_df)

        # Add locus tag and contig details to the RGI dataframes
        rgi_dataframes[rgi_filename] = make_RGI_df_contig_col(rgi_df)
        rgi_dataframes[rgi_filename]['Locus_Tag'], rgi_dataframes[rgi_filename]['Best_Hit_ARO'] = manipulate_GBK_contigs(rgi_df, rgi_filename)

        # Preprocess all bracketed columns to remove brackets
        for col in rgi_df:
            if type(rgi_df[col][0]) is list or isinstance(rgi_df[col][0], str):
                rgi_df[col] = rgi_df[col].apply(lambda x: strip_brackets(x))

    # 3) Divide GBK file data into contigs using group_by_contig to ensure we find neighbours on the same contig
    contig_dict = {}
    gbk_filenames = [get_filename(gbk_filepath) for gbk_filepath in gbk_filepaths]
    for gbk in gbk_filenames:
        for genome_id, gbk_df in gbk_dataframes.items():
            contig_dict[genome_id] = group_by_contig(gbk_df)

    # 4) Get unique ARO terms and union of all unique ARO terms
    ARO_best_terms, ARO_union = find_union_AMR_genes(rgi_dataframes)

    # 5) Get AMR dict for all genomes
    genomes_AMR_dict = {}
    for AMR_gene in ARO_union:
        # Make entry for gene occurrences
        genomes_AMR_dict[AMR_gene] = make_AMR_dict(rgi_dataframes, AMR_gene)

    # 6) Extract AMR neighborhoods and store them in neighborhood dataframes
    neighborhoods, flags = get_all_AMR_gene_neighborhoods(genomes_AMR_dict, gbk_dataframes, ARO_union, num_neighbours)

    print("Num neighborhoods: {}".format(len(neighborhoods)))
    print(neighborhoods)

    # 7) Get the locus, protein, and gene name details for each neighborhood respectively for FNA file creation
    locuses, protein_seqs, gene_names = get_neighborhood_data(neighborhoods)

    # 8) Save neighborhoods to FNA files needed for All-vs-all BLAST results later
    for AMR_gene, AMR_neighborhood_dict in neighborhoods.items():
        for genome_id, neighborhood in AMR_neighborhood_dict.items():
            # Make output subdirectory for the gene
            output_path = 'output/' + AMR_gene
            check_output_path(output_path)

            # Make file with genome ID as filename
            write_AMR_neighborhood_to_FNA(neighborhoods, AMR_gene, genome_id, locuses, protein_seqs, output_path)

    execution_time = (time.time() - start_time)
    print("Execution time in seconds: {}".format(str(execution_time)))