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
from utils import get_filename, check_output_path, strip_brackets, generate_alphanumeric_string
from filtering import filter_neighborhoods, write_filtered_genomes_textfile


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
    rgi_file_names = set(os.path.basename(file).strip('_rgi.txt') for file in rgi_filepaths)
    gbk_file_names = set(os.path.basename(file).strip('.gbk') for file in gbk_filepaths)

    #assert rgi_file_names == gbk_file_names, "Error: mismatch occurred between RGI and GBK file names."

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


def adjust_RGI_df_orientation(rgi_df):
    """
    Replaces RGI orientation symbols (-, +) with GBK convention (-1, +1) for easier cross-comparison
    """
    rgi_df['Orientation'] = rgi_df['Orientation'].map({'-': '-1',
                                                       '+': '+1'},
                                                      na_action=None)
    return rgi_df


def swap_neighborhood_orientation(df):
    """
    Reverses neighborhood representation. Used to ensure all AMR genes have same orientation in downstream gene order
    visualizations.
    """
    df['Gene_Strand'] = df['Gene_Strand'].map({-1: +1,
                                               +1: -1},
                                               na_action=None)

    return df

def clean_AMR_gene_name(temp_name):
    """
    Shortens AMR gene name and removes restricted characters that can cause issues in downstream analysis.
    """
    if len(temp_name) > 1:
        name = temp_name[0] + '_' + temp_name[1]
        # name = temp_name[0][0] + temp_name[1][0] + "_" + temp_name[2]
    else:
        name = temp_name[0]

    # Remove characters
    AMR_gene_name = name.replace('(', '').replace(')', '').replace('/', '_').replace('-','_').strip('_')

    return AMR_gene_name


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

        name = clean_AMR_gene_name(temp_name)
        RGI_names.append(name)

        locus_tags.append(genome_name + '(' + name + ')' + dataframe['Cut_Off'][index][0] +
                          '_' + str(dataframe['Best_Identities'][index]))

    return locus_tags, RGI_names


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
        print(genome_id)

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
                neighbor = contig_df.iloc[downstream_indices[i]]
                downstream_neighbors = downstream_neighbors.append(neighbor)
            except IndexError:
                print("Contig end found at position -{} downstream.".format(i + 1))
                neighborhood_indices.append(i + 1)

        # If there was no contig end, append default N size to neighborhood_indices
        if len(neighborhood_indices) == 0:
            neighborhood_indices.append(-neighborhood_size)

        # Get upstream neighbors
        upstream_indices = [gene_index + index for index in range(1, neighborhood_size + 1)]
        upstream_neighbors = pd.DataFrame(columns=['Gene_Start', 'Gene_End', 'Gene_Strand', 'Locus_Tag', 'Gene_Name',
                                                   'Product', 'Protein_Sequence', 'Contig_Name'])
        for i in range(len(upstream_indices)):
            contig_end_found = False
            try:
                neighbor = contig_df.iloc[upstream_indices[i]]
                upstream_neighbors = upstream_neighbors.append(neighbor)
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

        # Keep track of start and stop indices for each neighborhood
        contig_end_flags = {}

        for genome, data in AMR_dict.items():

            # Get start index of the focal AMR gene from the RGI dataframe
            start_vals = AMR_dict[genome]['Start']
            start_vals_list = list(start_vals)
            start_index = start_vals_list[0]

            # Make gene neighborhood dataframe for each genome for the focal gene, AMR_gene
            try:
                neighborhood_df, indices = make_AMR_gene_neighborhood_df(GBK_df_dict, genome, start_index, AMR_gene,
                                                                         neighborhood_size)
                neighborhood_df.reset_index(drop=True, inplace=True)

                if len(neighborhood_df) > 1:
                    neighbors[genome] = neighborhood_df
                    contig_end_flags[genome] = indices
            except TypeError:
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
                output_fasta.write('>' + genome_id + '_{}\n'.format(locus_output))
                output_fasta.write('{}\n'.format(protein_output))
        output_fasta.close()


def delete_low_occurring_genes(AMR_gene_dict, num_genomes, cutoff_percentage=0.25):
    """
    Removes keys of AMR genes not present in cutoff percentage of genomes.
    - AMR_gene_dict should be the dictionary of one type of AMR gene, with entries representing
    its occurrences within genomes.
    - Cutoff percentage is a float between 0 and 1 (e.g., 0.25 means only genes present in min 35% of genomes are kept).
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
        perfect_count = len(rgi_df[rgi_df["Cut_Off"] == "Perfect"])
        strict_count = len(rgi_df[rgi_df["Cut_Off"] == "Strict"])
        loose_count = len(rgi_df[rgi_df["Cut_Off"] == "Loose"])
        amr_gene_statistics[genome] = [perfect_count, strict_count, loose_count]

    return amr_gene_statistics


def write_summary_file(output_dir, gbk_dataframes, neighborhood_size, ARO_union, amr_gene_statistics):
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
        outfile.write('Number of AMR gene models extracted: {}'.format(len(ARO_union)) + '\n\n')
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        outfile.write('AMR Gene Details Per Genome' + '\n')
        outfile.write('----------------------------------------------------------------------------------------' + '\n')
        outfile.write('GENOME_ID' + '\t' + 'PERFECT_HITS' + '\t' + 'STRICT_HITS' + '\t' + 'LOOSE_HITS' + '\n')
        for genome, counts in amr_gene_statistics.items():
            outfile.write('{genome}\t{perfect}\t{strict}\t{loose}\n'.format(genome=get_filename(genome),
                                                                            perfect=counts[0],
                                                                            strict=counts[1],
                                                                            loose=counts[2]))


def reverse_start_end(df, n_start, n_stop):
    """
    Assuming other aspects of gene neighborhood df have been reversed (e.g. orientation, order, downstream/upstream),
    modifies each neighborhood gene's start and stop index to reflect these changes.
    """
    for gene in df.index:
        start = df['Gene_Start'][gene]
        end = df['Gene_End'][gene]
        gene_length = end - start

        # Subtract end coord
        adjusted_start = start - n_stop
        adjusted_end = end - n_stop

        # Flip values and swap
        new_start = abs(adjusted_end)
        new_end = abs(adjusted_start)

        # Update indices
        df.loc[gene, ['Gene_Start']] = [new_start + n_start]
        df.loc[gene, ['Gene_End']] = [new_end + n_start]

    return df


def reverse_df(df, n_start, n_stop, amr_gene_index):
    """
    Given a neighborhood dataframe we want to reverse, performs the following:
    (i) Reverses orientation of each gene.
    (ii) Reverses order genes appear in (e.g. downstream genes [A, B, C] -> [C, B, A]).
    (iii) Reverses upstream/downstream genes (e.g. [upstream, AMR, downstream] -> [downstream, AMR, upstream]).
    (iv) Modify gene start and end positions to reflect neighborhood reversal.
    """
    # (iv) Reverse gene start/end positions to reflect reversal
    neighborhood_df = reverse_start_end(df, n_start, n_stop)

    # (i) Swap gene orientations
    swapped_df = swap_neighborhood_orientation(neighborhood_df)
    swapped_df.reset_index(drop=True, inplace=True)
    sorted_df = swapped_df.sort_values(by='Gene_Start')
    sorted_df.reset_index(drop=True, inplace=True)

    return swapped_df.copy(deep=True)


def make_neighborhood_JSON_data(AMR_gene_neighborhoods_dict, AMR_gene):
    """
    Creates dictionary of data required to write neighborhoods JSON file.
    This function should be run on the complete set of neighborhoods (i.e. for all AMR genes).
    """
    # Keeps track of data for clusters, links, and groups
    neighborhood_json_data = {}

    # CLUSTER DATA
    cluster_data = {}
    unique_genes = []

    # Genome gene contig dict
    genome_contigs = {}

    for genome_id, neighborhood_df in AMR_gene_neighborhoods_dict.items():

        neighborhood_data = {}
        contigs_dict = {}

        # cluster_uids: Assign a random alphanumeric string of 36 characters for uid (as per clustermap.js examples),
        neighborhood_data["uid"] = generate_alphanumeric_string(36)
        neighborhood_data["name"] = genome_id

        loci_data = {}
        loci_data["uid"] = generate_alphanumeric_string(36)
        loci_data["name"] = genome_id

        # Retain neighborhood start/end coords
        df = neighborhood_df.copy(deep=True)
        loci_data["start"] = df['Gene_Start'].min()
        loci_data["end"] = df['Gene_End'].max()

        # cluster_loci: Key: uid (gene index in neighborhood), Values: (arr) gene name, start coord, end coord, strand
        genes_dict = {}

        indices = neighborhood_df.index.values.tolist()

        # Reverse the neighborhood representation if the focal gene is not oriented correctly
        amr_gene_index = int((len(indices) - 1) / 2)

        if neighborhood_df['Gene_Strand'][amr_gene_index] == -1:
            neighborhood_df = reverse_df(neighborhood_df, loci_data["start"], loci_data["end"], amr_gene_index)

        for i in neighborhood_df.index:
            gene_data = {}
            gene_data["uid"] = loci_data["uid"] + '-' + str(i)
            gene_data["name"] = neighborhood_df['Gene_Name'][i]
            if neighborhood_df['Gene_Name'][i].replace("'", "") not in unique_genes:
                unique_genes.append(neighborhood_df['Gene_Name'][i].replace("'", ""))
            gene_data["start"] = neighborhood_df['Gene_Start'][i]
            gene_data["end"] = neighborhood_df['Gene_End'][i]
            gene_data["strand"] = neighborhood_df['Gene_Strand'][i]

            # Keep track of contig data
            contigs_dict[loci_data["uid"] + '-' + str(i)] = neighborhood_df['Locus_Tag'][i]

            genes_dict[i] = gene_data

        loci_data["genes"] = genes_dict
        neighborhood_data["loci"] = loci_data
        cluster_data[genome_id] = neighborhood_data
        genome_contigs[genome_id] = contigs_dict

    neighborhood_json_data['clusters'] = cluster_data

    # LINK DATA
    links_data = {}

    # GROUPS DATA
    groups_data = {}
    ind = 0
    for i in range(len(unique_genes)):

        gene = unique_genes[i]
        genomes = list(cluster_data.keys())

        # Ignore unidentified genes
        if gene != 'UID':

            group_data = {}
            gene_group_list = []

            for genome_1, genome_2 in itertools.combinations(genomes, 2):

                link_data = {}

                genome_1_presence = False
                genome_2_presence = False
                genome_1_key = 0
                genome_2_key = 0

                for id in cluster_data[genome_1]["loci"]["genes"].keys():
                    if gene in cluster_data[genome_1]["loci"]["genes"][id]["name"]:
                        genome_1_presence = True
                        genome_1_key = id
                        break

                for id in cluster_data[genome_2]["loci"]["genes"].keys():
                    if gene in cluster_data[genome_2]["loci"]["genes"][id]["name"]:
                        genome_2_presence = True
                        genome_2_key = id
                        break

                if genome_1_presence and genome_2_presence:

                    # Add link
                    target = {}
                    target["uid"] = cluster_data[genome_1]["loci"]["genes"][genome_1_key]["uid"]
                    target["name"] = cluster_data[genome_1]["loci"]["genes"][genome_1_key]["name"].replace("'", "")

                    query = {}
                    query["uid"] = cluster_data[genome_2]["loci"]["genes"][genome_2_key]["uid"]
                    query["name"] = cluster_data[genome_2]["loci"]["genes"][genome_2_key]["name"].replace("'", "")

                    # Add to group
                    if cluster_data[genome_1]["loci"]["genes"][genome_1_key]["uid"] not in gene_group_list:
                        gene_group_list.append(cluster_data[genome_1]["loci"]["genes"][genome_1_key]["uid"])
                    if cluster_data[genome_2]["loci"]["genes"][genome_2_key]["uid"] not in gene_group_list:
                        gene_group_list.append(cluster_data[genome_2]["loci"]["genes"][genome_2_key]["uid"])

                    contig_1 = genome_contigs[genome_1][target["uid"]].replace("'", "")
                    contig_2 = genome_contigs[genome_2][query["uid"]].replace("'", "")
                    link_data["uid"] = str(genome_1) + '_' + contig_1 + '-' + str(genome_2) + '_' + contig_2
                    link_data["target"] = target
                    link_data["query"] = query

                    # TO DO: replace with BLAST P.I
                    link_data["identity"] = 0.5

                    links_data[ind] = link_data
                    ind += 1

            group_data["uid"] = 'group' + str(i)
            group_data["label"] = gene
            group_data["genes"] = gene_group_list
            groups_data[i] = group_data

    neighborhood_json_data['links'] = links_data
    neighborhood_json_data['groups'] = groups_data

    return neighborhood_json_data


def write_neighborhood_JSON(neighborhood_json_data, AMR_gene, output_path, surrogates=False):
    """
    Creates JSON file containing neighborhood data for an AMR gene for the genomes under analysis of the following form
    (where the first N neighbor genes represent the N neighbors downstream of the target gene, and the last N
    neighbor genes represent the N neighbors upstream of the target gene).
    """
    # Write JSON format
    out_path = output_path + '/JSON'
    check_output_path(out_path)

    try:
        cluster_data = neighborhood_json_data['clusters']
        links_data = neighborhood_json_data['links']
        groups_data = neighborhood_json_data['groups']
    except KeyError:
        print("Neighborhood data for {g} was not found.".format(g=AMR_gene))
        return

    if surrogates:
        output_file_path = out_path + '/' + AMR_gene + '_surrogates.json'
    else:
        output_file_path = out_path + '/' + AMR_gene + '.json'

    with open(output_file_path, 'w') as outfile:
        outfile.write('{\n')

        # ----------------------------------------------- CLUSTERS -----------------------------------------------------
        outfile.write('\t"clusters": [\n')

        cluster_index = 0
        num_clusters = len(cluster_data.keys())
        for genome_id in cluster_data.keys():

            outfile.write('\t\t{\n')

            # Neighborhood uid and name
            outfile.write('\t\t\t"uid": ' + '"' + cluster_data[genome_id]["uid"] + '",\n')
            outfile.write('\t\t\t"name": ' + '"' + genome_id + '",\n')

            # Neighborhood details: uid, name, start, stop, end, and gene data
            outfile.write('\t\t\t"loci": [\n')
            outfile.write('\t\t\t\t{\n')
            outfile.write('\t\t\t\t\t"uid": ' + '"' + cluster_data[genome_id]["loci"]["uid"] + '",\n')
            outfile.write('\t\t\t\t\t"name": ' + '"' + cluster_data[genome_id]["loci"]["name"] + '",\n')
            outfile.write('\t\t\t\t\t"start": ' + str(cluster_data[genome_id]["loci"]["start"]) + ',\n')
            outfile.write('\t\t\t\t\t"end": ' + str(cluster_data[genome_id]["loci"]["end"]) + ',\n')

            # Neighborhood genes data
            outfile.write('\t\t\t\t\t"genes": [\n')

            index = 0
            neighborhood_length = len(cluster_data[genome_id]["loci"]["genes"])
            for gene_key, gene_data in cluster_data[genome_id]["loci"]["genes"].items():

                # Opening bracket
                outfile.write('\t\t\t\t\t\t{\n')

                # Gene data
                outfile.write(
                    '\t\t\t\t\t\t\t"uid": ' + '"' + cluster_data[genome_id]["loci"]["genes"][gene_key]["uid"] + '",\n')
                outfile.write(
                    '\t\t\t\t\t\t\t"name": ' + '"' + cluster_data[genome_id]["loci"]["genes"][gene_key]["name"].replace(
                        "'", "") + '",\n')
                outfile.write('\t\t\t\t\t\t\t"start": ' + str(
                    cluster_data[genome_id]["loci"]["genes"][gene_key]["start"]) + ',\n')
                outfile.write(
                    '\t\t\t\t\t\t\t"end": ' + str(cluster_data[genome_id]["loci"]["genes"][gene_key]["end"]) + ',\n')
                outfile.write('\t\t\t\t\t\t\t"strand": ' + str(
                    cluster_data[genome_id]["loci"]["genes"][gene_key]["strand"]) + '\n')

                # Closing bracket
                if index != len(cluster_data[genome_id]["loci"]["genes"]) - 1:
                    outfile.write('\t\t\t\t\t\t},\n')
                    index += 1
                else:
                    outfile.write('\t\t\t\t\t\t}\n')

            # Closing brackets for clusters
            outfile.write('\t\t\t\t\t]\n')
            outfile.write('\t\t\t\t}\n')
            outfile.write('\t\t\t]\n')
            if cluster_index != num_clusters - 1:
                outfile.write('\t\t},\n')
                cluster_index += 1
            else:
                outfile.write('\t\t}\n')

        # Clusters data final closing bracket
        outfile.write('\t],\n')

        # ----------------------------------------------- LINKS --------------------------------------------------------
        outfile.write('\t"links": [\n')

        link_index = 0
        num_links = len(links_data.keys())
        last_link_flag = False
        for link in links_data.keys():

            try:
                uid = links_data[link]["target"]["uid"]
                outfile.write('\t\t{\n')
                outfile.write('\t\t\t"uid": "' + links_data[link]["uid"] + '",\n')
                outfile.write('\t\t\t"target": {\n')
                outfile.write('\t\t\t\t"uid": "' + links_data[link]["target"]["uid"] + '",\n')
                outfile.write('\t\t\t\t"name": "' + links_data[link]["target"]["name"] + '-1"\n')
                outfile.write('\t\t\t},\n')
                outfile.write('\t\t\t"query": {\n')
                outfile.write('\t\t\t\t"uid": "' + links_data[link]["query"]["uid"] + '",\n')
                outfile.write('\t\t\t\t"name": "' + links_data[link]["query"]["name"] + '-2"\n')
                outfile.write('\t\t\t},\n')
                outfile.write('\t\t\t"identity": "' + str(links_data[link]["identity"]) + '"\n')

            except KeyError:
                pass

            if link_index != num_links - 1:
                outfile.write('\t\t},\n')

            else:
                outfile.write('\t\t}\n')
                last_link_flag = True

            link_index += 1

            if last_link_flag:
                break

        # Links data final closing bracket
        outfile.write('\t],\n')

        # ------------------------------------------------ GROUPS ------------------------------------------------------
        outfile.write('\t"groups": [\n')

        group_index = 0
        num_groups = len(groups_data.keys())
        for group in groups_data.keys():
            outfile.write('\t\t{\n')

            outfile.write('\t\t\t"uid": "' + groups_data[group]["uid"] + '",\n')
            outfile.write('\t\t\t"label": "' + groups_data[group]["label"] + '",\n')
            outfile.write('\t\t\t"genes": ' + json.dumps(groups_data[group]["genes"]) + '\n')

            if group_index != num_groups - 1:
                outfile.write('\t\t},\n')
                group_index += 1
            else:
                outfile.write('\t\t}\n')

        outfile.write('\t]\n')

        # Final closing bracket
        outfile.write('}\n')

        return


def make_AMR_gene_HTML(AMR_genes_list, sample_data_path, out_path):
    """
    For each AMR gene for which a JSON file was created, generates an accompanying HTML file for rendering its gene
    order visualization using clustermap with. This is done for each gene individually to
    """

    for AMR_gene in AMR_genes_list:

        # Make HTML index file with appropriate JSON
        with open(out_path + '/' + AMR_gene + '.html', 'w') as html_outfile, open(
                sample_data_path + '/index.html') as template:
            for line in template:
                html_outfile.write(line)

            html_outfile.write('\n')
            html_outfile.write('\t\td3.json("' + AMR_gene + '.json"' + ')\n')
            html_outfile.write('\t\t\t.then(data => {\n')
            html_outfile.write('\t\t\t\tdiv.selectAll("div")\n')
            html_outfile.write('\t\t\t\t\t.data([data])\n')
            html_outfile.write('\t\t\t\t\t.join("div")\n')
            html_outfile.write('\t\t\t\t\t.call(chart)\n\n')
            html_outfile.write('\t\t\t\tlet svg = div.select("svg")\n')
            html_outfile.write('\t\t\t\td3.select("#btn-save-svg")\n')
            html_outfile.write('\t\t\t\t\t.on("click", () => {\n')
            html_outfile.write('\t\t\t\t\t\tconst blob = serialise(svg)\n')
            html_outfile.write('\t\t\t\t\t\tdownload(blob, "clinker.svg")\n')
            html_outfile.write('\t\t\t\t\t})\n')
            html_outfile.write('\t\t\t})\n')
            html_outfile.write('\t</script>\n')
            html_outfile.write('</html>')

        # Make surrogate HTML
        with open(out_path + '/' + AMR_gene + '_surrogates.html', 'w') as html_outfile, open(
                sample_data_path + '/index.html') as template:
            for line in template:
                html_outfile.write(line)

            html_outfile.write('\n')
            html_outfile.write('\t\td3.json("' + AMR_gene + '_surrogates.json"' + ')\n')
            html_outfile.write('\t\t\t.then(data => {\n')
            html_outfile.write('\t\t\t\tdiv.selectAll("div")\n')
            html_outfile.write('\t\t\t\t\t.data([data])\n')
            html_outfile.write('\t\t\t\t\t.join("div")\n')
            html_outfile.write('\t\t\t\t\t.call(chart)\n\n')
            html_outfile.write('\t\t\t\tlet svg = div.select("svg")\n')
            html_outfile.write('\t\t\t\td3.select("#btn-save-svg")\n')
            html_outfile.write('\t\t\t\t\t.on("click", () => {\n')
            html_outfile.write('\t\t\t\t\t\tconst blob = serialise(svg)\n')
            html_outfile.write('\t\t\t\t\t\tdownload(blob, "clinker.svg")\n')
            html_outfile.write('\t\t\t\t\t})\n')
            html_outfile.write('\t\t\t})\n')
            html_outfile.write('\t</script>\n')
            html_outfile.write('</html>')


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

    # Make output directory if non-existent
    check_output_path(output_path)

    print("Processing GBK and RGI files...")

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

    print(gbk_dataframes)

    # 2) Get RGI dataframes and do all required preprocessing for later manipulation and comparison
    # Contains a dataframe for each RGI file, where keys are filenames
    rgi_dataframes = {}

    for rgi_filepath in rgi_filepaths:
        rgi_filename = get_filename(rgi_filepath, rgi=True)
        rgi_df = make_RGI_dataframe(rgi_filepath)
        rgi_dataframes[rgi_filename] = rgi_df

    amr_statistics = get_AMR_gene_statistics(rgi_dataframes)

    for rgi_filename, rgi_df in rgi_dataframes.items():

        # Replace each RGI dataframe with one with the proper DNA orientation representation to match GBKs
        rgi_dataframes[rgi_filename] = adjust_RGI_df_orientation(rgi_df)

        # Add locus tag and contig details to the RGI dataframes
        rgi_dataframes[rgi_filename] = make_RGI_df_contig_col(rgi_df)
        rgi_dataframes[rgi_filename]['Locus_Tag'], rgi_dataframes[rgi_filename][
            'Best_Hit_ARO'] = manipulate_GBK_contigs(rgi_df, rgi_filename)

        # Preprocess all bracketed columns to remove brackets
        for col in rgi_df:
            if type(rgi_df[col][0]) is list or isinstance(rgi_df[col][0], str):
                try:
                    rgi_df[col] = rgi_df[col].apply(lambda x: strip_brackets(x))
                except AttributeError:
                    pass

    # 3) Divide GBK file data into contigs using group_by_contig to ensure we find neighbors on the same contig
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

    genomes_AMR_dict_filtered = delete_low_occurring_genes(genomes_AMR_dict, len(gbk_filepaths), cutoff_percent)

    print("Extracting gene neighborhood data for neighborhoods of size {}...".format(num_neighbors))

    # 6) Extract AMR neighborhoods and store them in neighborhood dataframes
    print(genomes_AMR_dict_filtered.keys())
    neighborhoods, neighborhood_indices = get_all_AMR_gene_neighborhoods(genomes_AMR_dict_filtered, gbk_dataframes,
                                                                         ARO_union, num_neighbors)

    # 7) Get the locus, protein, and gene name details for each neighborhood respectively for FNA file creation
    locuses, protein_seqs, gene_names = get_neighborhood_data(neighborhoods, num_neighbors)

    # 8) Save neighborhoods to FNA files needed for All-vs-all BLAST results later
    print("Generating neighborhood FNA files...")
    for AMR_gene, genome_neighborhoods in neighborhoods.items():

        # Make output subdirectory for the gene
        gene_name = AMR_gene.replace('(', '_').replace(')', '_').replace('-', '').replace("'", "")
        out_path = output_path + '/fasta/' + gene_name
        check_output_path(out_path)

        for genome_id in genome_neighborhoods.keys():
            # Make file with genome ID as filename
            write_AMR_neighborhood_to_FNA(genome_neighborhoods, AMR_gene, locuses, protein_seqs, out_path)

    # 9) Make neighborhoods summary textfile in output dir
    print("Making extraction summary file...")
    write_summary_file(output_path, gbk_dataframes, num_neighbors, ARO_union, amr_statistics)

    # 10) Save gene neighborhoods and indices in textfile: needed for JSON representations downstream
    for AMR_gene, neighborhood_data in neighborhoods.items():

        # Filter neighborhoods
        filtered_neighborhoods = filter_neighborhoods(neighborhood_data)
        surrogates = list(filtered_neighborhoods.keys())
        filtered_neighborhoods_dict = {key : value for (key, value) in neighborhood_data.items() if key in surrogates}
        write_filtered_genomes_textfile(filtered_neighborhoods, AMR_gene, output_path)

        # Get data needed to write JSON files
        neighborhood_JSON_dict = make_neighborhood_JSON_data(neighborhood_data, AMR_gene)
        filtered_neighborhood_JSON_dict = make_neighborhood_JSON_data(filtered_neighborhoods_dict, AMR_gene)

        # Create JSON file
        write_neighborhood_JSON(neighborhood_JSON_dict, AMR_gene, output_path)
        write_neighborhood_JSON(filtered_neighborhood_JSON_dict, AMR_gene, output_path, True)

    sample_data_path = '../sample_data'
    make_AMR_gene_HTML(neighborhoods.keys(), sample_data_path, output_path)

    neighborhoods_dict = {}
    for AMR_gene, genome in neighborhoods.items():
        genome_dict = {}
        for genome, neighborhood_df in genome.items():
            neighborhood_dict = neighborhood_df.to_dict()
            genome_dict[genome] = neighborhood_dict
        neighborhoods_dict[AMR_gene] = genome_dict

    with open(output_path + '/' + 'neighborhood_data.txt', 'w') as outfile:
        for AMR_gene, genome in neighborhoods_dict.items():
            print(AMR_gene, file=outfile)
            for genome, neighborhood_data in genome.items():
                print('\t{}'.format(genome), file=outfile)
                for col, row_data in neighborhood_data.items():
                    print('\t\t{}:'.format(col), file=outfile)
                    for row_id, row_val in row_data.items():
                        print('\t\t\t{}: {}'.format(row_id, row_val), file=outfile)
                    print(file=outfile)
                print(file=outfile)
                print("---GENOME-DATA-END---", file=outfile)
            print("---GENE-DATA-END---", file=outfile)
            print(file=outfile)

    with open(output_path + '/' + 'neighborhood_indices.txt', 'w') as outfile:
        outfile.write(str(neighborhood_indices))

    print("Neighborhood extraction complete.")
    check_output_path(output_path + '/blast')


def main(args=None):
    args = parse_args(args)
    extract_neighborhoods(args.RGI_PATH, args.GBK_PATH, args.OUTPUT_PATH, args.n, args.p)


if __name__ == '__main__':
    sys.exit(main())
