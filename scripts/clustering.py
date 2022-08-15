#!/usr/bin/env python

"""
Functions to cluster neighborhoods identified using get_neighborhoods.py.
"""

import numpy as np
import pandas as pd
import itertools
import glob
import os
import shutil
from scipy.cluster.hierarchy import average, fcluster
from scipy import sparse
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN
# import markov_clustering as mcl
from extraction import get_filename, check_output_path
from utils import get_full_filepaths, remove_files
from scoring import get_normalized_bitscores

def run_BLAST(neighborhoods_filepaths, root_path, output_path):
    """
    Given filepaths for all genome neighborhoods of an AMR gene X, runs All-vs-All BLAST
    and stores the results in output/BLAST/X.
    """
    # BLAST each neighborhood against every other neighborhood once
    for neighborhood_fasta_1, neighborhood_fasta_2 in itertools.combinations(neighborhoods_filepaths, 2):

        print("Combo: {a} and {b}".format(a=neighborhood_fasta_1, b=neighborhood_fasta_2))

        # Get the filepaths for the neighborhoods to BLAST this iteration
        fasta_1_filename = get_filename(neighborhood_fasta_1)
        fasta_2_filename = get_filename(neighborhood_fasta_2)
        blast_output_filename = fasta_1_filename + '_' + fasta_2_filename + ".blast.txt"

        # Make output dir path and switch directories
        blast_output_path = output_path + '/' + blast_output_filename

        # Make BLAST db for the first neighborhood
        os.system('makeblastdb -in ' + neighborhood_fasta_1 + ' -parse_seqids -dbtype prot')

        # Run BLAST All-vs-All with specified e-value cutoff
        os.system('blastp -query ' + neighborhood_fasta_2 + ' -db ' + neighborhood_fasta_1 + \
                  ' -task blastp -outfmt 6 -max_hsps 1 -out ' + blast_output_filename)

        # Move file to correct directory in at output_path
        curr_dir = os.getcwd()
        curr_blast_path = os.path.join(curr_dir, blast_output_filename)
        shutil.move(curr_blast_path, blast_output_path)

def make_BLAST_df_dict(blast_gene_dir):
    """
    Returns a dictionary of BLAST dataframes, where a dataframe exists for each BLAST file and the name is the key.
    """
    #blast_filenames = [os.path.basename(filepath) for filepath in blast_filepaths]

    blast_dataframes = {}
    if os.path.isdir(blast_gene_dir):
        for filename in os.listdir(blast_gene_dir):
            #with open(os.path.join(blast_gene_dir, filename)) as blast_file:
                #data = blast_file.read()
            file = os.path.join(blast_gene_dir, filename)
            blast_file_df = pd.read_csv(file, index_col=False, sep='\t', names=['query_id', 'sub_id', 'PI', 'aln_len',
                                                                                'n_of_mismatches', 'gap_openings',
                                                                                'q_start', 'q_end', 's_start', 's_end',
                                                                                'Evalue', 'bitscore'])
            # Ensure correct datatypes for cols
            for col in ['PI', 'aln_len', 'n_of_mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end',
                        'Evalue', 'bitscore']:
                blast_file_df[col] = blast_file_df[col].astype(int)
            blast_dataframes[filename] = blast_file_df

    return blast_dataframes

def filter_BLAST_results(blast_gene_dict, cutoff=70):
    """
    Retains only protein sequences that share more than 70% overall sequence identity (i.e., homologs) in BLAST dict.
    """
    for blast_file, blast_df in blast_gene_dict.items():
        # Filter rows to drop all rows with percent identity lower than specified threshold
        df = blast_df[blast_df['PI'] < cutoff]
        blast_gene_dict[blast_file] = df

    return blast_gene_dict

def remove_BLAST_duplicates(blast_gene_dict):
    """
    Given a BLAST file dataframe, removes duplicate entries for same gene by keeping only the highest bitscore.
    """
    for blast_file, blast_df in blast_gene_dict.items():
        blast_df = blast_df.sort_values('bitscore', ascending=False)
        blast_df = blast_df.drop_duplicates(['query_id'])
        blast_df.reset_index(drop=True, inplace=True)

    return blast_gene_dict

def get_neighborhood_similarity_matrix(neighborhoods_dict, num_neighbours):
    """
    Given a neighborhood dictionary for a given drug obtained from get genes,
    calculate similarity scores and return a similarity matrix. (WIP!)
    """
    similarity_matrix_dict = {}
    for key, item in neighborhoods_dict.items():
        similarity_matrix = []

        for genome in neighborhoods_dict[key].keys():
            gene = neighborhoods_dict[key][genome]

            for genome_to_compare, i in item.items():
                temp_data = []
                sum_1 = 0
                sum_2 = 0

                gene.reset_index(drop=True, inplace=True)

                for index in range(len(gene)):
                    if gene['sub_id'][index] in i:
                        temp_data.append((gene['query_data'][index],
                                          gene['sub_id'][index],
                                          gene['normalized_bitscore'][index]))

                temp_df = remove_BLAST_duplicates(temp_data)
                sum_1 = round(temp_df['bitscore'].sum(), 3)

                # SINGLE CONTIG END FLAG DICT HERE (block to handle factoring in contig ends)
                #if mult_contig_end_dict[key][genome] == 1:
                #    len_neighborhood = len(neighborhoods_dict[key][genome])
                #    len_diff = len_neighborhood - len(temp_df)
                #    sum_2 = sum_1 + (len_neighborhood - len_diff - len(temp_df))

                #elif mult_contig_end_dict[key][genome_to_compare] == 1:
                #    len_neighborhood = len(neighborhoods_dict[key][genome_to_compare])
                #    len_diff = len_neighborhood - len(temp_df)
                #    sum_2 = sum_1 + (len_neighborhood - len_diff - len(temp_df))

                #else:
                #    sum_2 = sum_1
                sum_2 = sum_1
                similarity_matrix.append((genome, genome_to_compare, round(sum_1, 3)))

        similarity_matrix_dict[key] = similarity_matrix

    return similarity_matrix_dict

def sim_matrix_to_symmetric_matrix(similarity_matrix):
    """
    Uses the averages of scores from similarity matrix to convert it into a symmetric matrix for clustering.
    """
    symmetric_matrix = []
    for row in range(len(similarity_matrix)):
        for col in range(len(similarity_matrix)):
            if similarity_matrix[col][1] == similarity_matrix[row][0] and \
                    similarity_matrix[col][0] == similarity_matrix[row][1]:
                # Get rounded score approximated to 3 decimal points
                average_score = round((similarity_matrix[row][2] + similarity_matrix[col][2]) / 2, 3)

        symmetric_matrix.append((similarity_matrix[row][0], similarity_matrix[row][1], average_score))

    return symmetric_matrix

def get_sparse_matrix(np_matrix):
    """
    Converts a numpy distance matrix to a scipy sparse matrix (important for use with MCL clustering).
    """
    sparse_matrix = sparse.csr_matrix(np_matrix)

    return sparse_matrix

def UPGMA_clustering(np_distance_matrix):
    """
    Applies UPGMA clustering to neighborhood similarity or symmetric distance matrix.
    """
    # UPGMA clustering to obtain dendogram
    upgma_dendogram = average(np_distance_matrix)

    return upgma_dendogram

def agglomerative_clustering(np_distance_matrix):
    """
    Applies agglomerative clustering to a neighborhood similarity or symmetric distance matrix.
    """
    acl = AgglomerativeClustering().fit(np_distance_matrix)

    return acl, acl.labels_

def MCL_clustering(sparse_distance_matrix):
    """
    Performs MCL clustering on a neighborhood sparse similarity or symmetric distance matrix.
    """
    # Need to explore with the inflation parameter.. Documentation suggests trying 1.4, 2, 4, and 6
    result = mcl.run_mcl(sparse_distance_matrix)
    clusters = mcl.get_clusters(result)

    return clusters

def DBSCAN_clustering(np_distance_matrix):
    """
    Applies DBSCAN clustering to a neighborhood similarity or symmetric distance matrix.
    """
    dbscan = DBSCAN(eps=0.5)  # Need some form of hyperparameter tuning... Optuna, since > gridsearch?

    return dbscan, dbscan.labels_

if __name__ == '__main__':

    try:
        fasta_path = sys.argv[1]
    except IndexError:
        print("Path to FASTA input files not correctly specified.")

    try:
        output_path = sys.argv[2]
    except IndexError:
        print("Output path not specified.")

    try:
        num_neighbors = sys.argv[3]
    except IndexError:
        # Default value
        num_neighbors = 10

    # Get BLAST results for each AMR gene neighborhood set
    fasta_output_filepath = 'output/fasta'
    fasta_folder_paths = get_full_filepaths(fasta_output_filepath)

    curr_dir = os.getcwd()
    BLAST_output_dir = os.path.join(curr_dir, 'output/blast')
    check_output_path(BLAST_output_dir)

    # For each AMR gene, get the neighborhood set for all genomes for a given AMR gene
    for AMR_gene_subdir in fasta_folder_paths:

        # Remove non-FASTA files
        remove_files(AMR_gene_subdir, '.fasta')

        # Make output dir for AMR gene within blast output
        AMR_gene_name = get_filename(AMR_gene_subdir)

        output_path = BLAST_output_dir + '/' + AMR_gene_name
        if os.path.exists(output_path):
            pass
        else:
            check_output_path(output_path)
            # Get all genome neighborhood fasta paths
            root_path = fasta_output_filepath + '/' + AMR_gene_name
            full_genome_paths = get_full_filepaths(root_path)

            # All vs All BLAST each genome pair against each other once
            run_BLAST(full_genome_paths, root_path, output_path)

    # Make BLAST dataframe dictionary
    blast_folder_paths = get_full_filepaths(BLAST_output_dir)

    BLAST_df_dict = {}
    for blast_amr_subdir in blast_folder_paths:

        # Get BLAST filenames
        blast_filepaths = get_full_filepaths(blast_amr_subdir)
        blast_filenames = [get_filename(path) for path in blast_filepaths]

        # For each gene, get dictionary of all blast results above %ID criterion
        blast_gene_dict = make_BLAST_df_dict(blast_amr_subdir)
        filtered_blast_dict = filter_BLAST_results(blast_gene_dict)

        # Remove duplicates
        blast_dict = remove_BLAST_duplicates(filtered_blast_dict)

        # Calculate normalized bitscores and save as column
        for filename, blast_df in blast_dict.items():
            normalized_bitscores = get_normalized_bitscores(blast_df)
            blast_df['normalized_bitscores'] = normalized_bitscores
            blast_dict[filename] = blast_df

        # Store data
        AMR_gene = get_filename(blast_amr_subdir)
        BLAST_df_dict[AMR_gene] = filtered_blast_dict

    # Calculate similarity matrix for each gene
    similarity_matrices_dict = get_neighborhood_similarity_matrix(BLAST_df_dict, 10)

    # Make output dir for clustering results
    check_output_path(output_path + '/clustering/UPGMA')

    # UPGMA clustering
    for AMR_gene, matrix_data in similarity_matrices_dict.items():
        check_output_path(output_path + '/' + AMR_gene)
        similarity_matrix = np.asmatrix(np.array(matrix_data))
        upgma_clusters = UPGMA_clustering(similarity_matrix)


