#!/usr/bin/env python

"""
Functions to cluster neighborhoods identified using extraction.py.
"""
import argparse
import itertools
import os
import shutil
import optuna
import json
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster import hierarchy
from scipy.spatial.distance import euclidean, pdist, squareform
from sklearn.cluster import DBSCAN
from sklearn.metrics import calinski_harabasz_score
from sklearn.preprocessing import StandardScaler

import markov_clustering as mcl

from DBCV import DBCV
from scoring import get_normalized_bitscores
from utils import get_filename, check_output_path, get_full_filepaths, remove_files, generate_alphanumeric_string, \
                  write_clustermap_JSON_HTML, remove_defunct_clustermap_data
from visualization import plot_similarity_histogram, plot_distance_histogram, \
                          graph_UPGMA_clusters, draw_mcl_graph, graph_DBSCAN_clusters, \
                          plotly_pcoa, plotly_dendrogram, plotly_mcl_network


def parse_args(args=None):
    Description = "Cluster extracted AMR gene neighborhoods to compare conservation characteristics across genomes."
    Epilog = "Example usage: python clustering.py <ASSEMBLY_PATH> <FASTA_PATH> <BLAST_PATH> <OUTPUT_PATH> \
              -n <NEIGHBORHOOD_SIZE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('ASSEMBLY_PATH', metavar='asm_path', type=str,
                        help='Path to FA assembly files.')
    parser.add_argument('FASTA_PATH', metavar='fasta_path', type=str,
                        help='Path to directory containing neighborhood FASTA files.')
    parser.add_argument('BLAST_PATH', metavar='blast_path', type=str,
                        help='Path to directory containing BLAST text files of BLASTed neighborhoods.')
    parser.add_argument('OUTPUT_PATH', metavar='output_path', type=str, help='Path to output directory where '
                                                                             'extracted neighborhood FASTA files will'
                                                                             ' be saved.')
    parser.add_argument('-n', metavar='-neighborhood_size', type=int, default=10, help='Neighborhood window size, i.e. '
                                                                             'number of genes considered upstream and '
                                                                             'downstream of focal gene.')
    parser.add_argument('-i', metavar='-mcl_inf', type=int, default=2, help='Inflation hyperparameter for Markov \
                                                                                                           clustering.')
    parser.add_argument('-e', metavar='-dbscan_eps', type=float, default=0.5, help='Inflation hyperparameter for \
                                                                                                    DBSCAN clustering.')
    parser.add_argument('-m', metavar='-dbscan_min', type=int, default=5, help='Minimum samples hyperparameter for \
                                                                                                    DBSCAN clustering.')

    return parser.parse_args(args)


def run_DIAMOND_BLAST(neighborhoods_filepaths, dir_output_path, root_path=None):
    """
    Given filepaths for all genome neighborhoods of an AMR gene X, runs All-vs-All BLAST
    and stores the results in output/blast/{AMR_gene
    """
    # BLAST each neighborhood against every other neighborhood once
    for fasta_1, fasta_2 in itertools.combinations(neighborhoods_filepaths, 2):
        print("Blasting All-vs-All for {a} vs {b}...".format(a=fasta_1, b=fasta_2))

        # Get the filepaths for the neighborhoods to BLAST this iteration
        fasta_1_filename = get_filename(fasta_1)
        fasta_2_filename = get_filename(fasta_2)

        blast_output_path = os.path.abspath(dir_output_path)

        blast_output_path = os.path.abspath(dir_output_path)

        # Blast each file against itself: need bitscore data for bitscore normalization later
        os.system('diamond makedb --in ' + fasta_1 + ' -d ' + fasta_1 + '.dmnd')
        os.system('diamond blastp -d ' + fasta_1 + '.dmnd --query ' + fasta_1 + ' --outfmt 6 --out ' \
                  + fasta_1_filename + '.dmnd.out')

        os.system('diamond makedb --in ' + fasta_2 + ' -d ' + fasta_2 + '.dmnd')
        os.system('diamond blastp -d ' + fasta_2 + '.dmnd --query ' + fasta_2 + ' --outfmt 6 --max-hsps 1 --out ' \
                  + fasta_2_filename + '.dmnd.out')

        # Make output dir path and switch directories
        final_output_path = os.path.abspath(dir_output_path)
        #final_output_path = os.path.abspath(dir_output_path) + '/' + blast_output_filename

        # Make BLAST db for the first neighborhood
        #os.system('diamond makedb --in ' + fasta_1 + ' -d ' + fasta_1 + '.dmnd')

        # Run BLAST All-vs-All with specified e-value cutoff
        blast_output_filename = fasta_1_filename + '_' + fasta_2_filename + ".dmnd.out"
        os.system('diamond blastp -d ' + fasta_1 + '.dmnd --query ' + fasta_2 + \
                  ' --outfmt 6 --max-hsps 1 --out ' + blast_output_filename)

        # Combine all three files into the BLAST All-vs-All of the two genomes against each other and delete originals
        read_files = [fasta_1_filename + '.dmnd.out', fasta_2_filename + '.dmnd.out']
        with open(blast_output_filename, 'a') as outfile:
            for file in read_files:
                with open(file, 'r') as infile:
                    outfile.write(infile.read())

                # Delete file after appending
                os.remove(file)

        # Move file to correct directory in at output_path
        curr_dir = os.getcwd()
        curr_blast_path = os.path.join(curr_dir, blast_output_filename)
        print("Moving {b} to {o}...".format(b=curr_blast_path, o=final_output_path))
        shutil.move(curr_blast_path, final_output_path)


def blast_neighborhoods(assembly_path, blast_path):
    """
    Calls DIAMOND BLAST on all FAA files.
    """
    # Get blast directory FAA file paths
    fa_file_paths = get_full_filepaths(assembly_path)
    run_DIAMOND_BLAST(fa_file_paths, blast_path)


def blast_neighborhoods(fasta_path, blast_path):
    """
    Calls DIAMOND BLAST on all gene neighborhood FASTA outputs.
    """
    # Get BLAST results for each AMR gene neighborhood set
    fasta_folder_paths = get_full_filepaths(fasta_path)
    # For each AMR gene, get the neighborhood set for all genomes for a given AMR gene
    for AMR_gene_subdir in fasta_folder_paths:
        # Remove non-FASTA files
        remove_files(AMR_gene_subdir, '.fasta')

        # Make output dir for AMR gene within blast output
        AMR_gene_name = get_filename(AMR_gene_subdir)
        dir_output_path = blast_path + '/' + AMR_gene_name

        # If BLAST files already exist for a given gene directory, skip it
        if os.path.exists(dir_output_path):
                pass
        else:
                check_output_path(dir_output_path)
                # Get all genome neighborhood fasta paths
                root_path = fasta_path + '/' + AMR_gene_name
                full_genome_paths = get_full_filepaths(root_path)

                # All vs All BLAST each genome pair against each other once
                run_DIAMOND_BLAST(full_genome_paths, root_path, dir_output_path)


def make_fasta_contig_dict(fasta_path, AMR_gene):
    """
    Given FASTA neighborhood file, creates a dictionary where keys correspond to indices from 0 to N,
    and values are contig ids.
    """
    fasta_dict = {}
    for fasta_file in os.listdir(fasta_path + '/' + AMR_gene):
        if fasta_file.endswith('.fasta'):
            with open(fasta_path + '/' + AMR_gene + '/' + fasta_file, 'r') as infile:
                data = infile.readlines()
                fasta_file_dict = {}
                index = 0
                for line in data:
                    if line.startswith('>'):
                        contig_id = line.strip().replace('>', '')
                        fasta_file_dict[index] = contig_id
                        index += 1
                genome_id = fasta_file.split('.fasta')[0]
                fasta_dict[genome_id] = fasta_file_dict

    return fasta_dict


def get_blast_dict_whole_genomes(fasta_path, blast_path):
    """
    Creates dictionary for all AMR genes where each key is an AMR gene name (str) and each value is the dictionary
    containing genome combinations that were blasted together as keys and their accompanying blast file data in a
    dataframe as values.
    """
    AMR_genes = os.listdir(fasta_path)

    BLAST_df_dict = {}
    for AMR_gene in AMR_genes:

        # Get the fasta file contents for each genome for that gene so contigs can be easily parsed
        fasta_dict = make_fasta_contig_dict(fasta_path, AMR_gene)

        # Make the BLAST dataframe using the contig dict
        AMR_dict = get_AMR_blast_df(AMR_gene, blast_path, fasta_path, fasta_dict)

        #AMR_dict = get_AMR_blast_dict(blast_path, AMR_gene)
        BLAST_df_dict[AMR_gene] = AMR_dict

    return BLAST_df_dict


def load_BLAST_file_df(filename):
    """
    Loads a BLAST file assuming all data columns for output format 6 are present in default order.
    Ensures columns are loaded with correct variable types.
    """
    # Initialize empty dataframe to hold desired contig data
    with open(filename, 'r') as infile:
        df = pd.read_csv(infile, header=None, sep='\t', names=['query_id', 'sub_id', 'PI', 'aln_len',
                                                               'n_of_mismatches', 'gap_openings',
                                                               'q_start', 'q_end', 's_start', 's_end',
                                                               'Evalue', 'bitscore'])
        # Ensure correct datatypes for cols
        for col in ['PI', 'aln_len', 'n_of_mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end',
                    'Evalue', 'bitscore']:
            df[col] = df[col].astype(int)

    return df


def get_contig_rows(contig, df, identical=False):
    """
    Given a contig and a dataframe containing contents of a BLAST file, returns a dataframe of all rows
    with that contig as either the query or sub.
    If identical, find rows where query and sub id are the same.
    """
    if identical:
        contig_rows = df.loc[((df['query_id'] == contig) & (df['sub_id'] == contig))]
    else:
        contig_rows = df.loc[((df['query_id'] == contig) | (df['sub_id'] == contig))]

    return contig_rows


def append_FASTA_neighborhood_contigs(filename, AMR_gene, contig_dict, identical=False):
    """
    Given a dictionary of contigs in a given neighborhood, where keys are int and values correspond to contigs,
    iterates through the contigs to obtain all relevant BLAST rows to a given dataframe.
    """
    # Dataframe to append to
    df = pd.DataFrame(columns=['query_id', 'sub_id', 'PI', 'aln_len', 'n_of_mismatches', 'gap_openings',
                               'q_start', 'q_end', 's_start', 's_end', 'Evalue', 'bitscore'])

    # Dataframe containing genome BLAST results
    data_df = load_BLAST_file_df(filename)

    if AMR_gene == 'RCP1' and filename == 'SAMEA1486355_SAMEA1466699.dmnd.out':
        data_df.to_csv('data_df.csv', index=True)

    for contig in contig_dict.keys():
        contig_rows = get_contig_rows(contig_dict[contig], data_df, identical=identical)
        df = df.append(contig_rows)
        df.reset_index(drop=True, inplace=True)

    del data_df

    return df


def get_AMR_blast_df(AMR_gene, blast_path, fasta_path, fasta_dict):
    """
    Loads whole genome BLAST into a dataframe with only relevant contig data rows.
    """
    # Identify all genomes the gene is present in using FASTA outputs
    present_genomes = [key.split('.')[0] for key in fasta_dict.keys()]

    AMR_dict = {}
    for genome_1, genome_2 in itertools.combinations(present_genomes, 2):

        # Initialize empty dataframe to hold desired contig data
        df = pd.DataFrame(columns=['query_id', 'sub_id', 'PI', 'aln_len', 'n_of_mismatches', 'gap_openings',
                                   'q_start', 'q_end', 's_start', 's_end', 'Evalue', 'bitscore'])

        # Load comparative BLAST file
        filename = ''
        genome_1_genome_2_file_name = blast_path + '/'  + genome_1 + '_' + genome_2 + '.dmnd.out'
        genome_2_genome_1_file_name = blast_path + '/' + genome_2 + '_' + genome_1 + '.dmnd.out'
        genome_1_file_name = blast_path + '/' + genome_1 + '_' + genome_1 + '.dmnd.out'
        genome_2_file_name = blast_path + '/' + genome_2 + '_' + genome_2 + '.dmnd.out'

        # Append all relevant BLAST rows for neighborhood contigs to df
        files = [genome_1_genome_2_file_name, genome_2_genome_1_file_name, genome_1_file_name, genome_2_file_name]
        contig_dicts = [fasta_dict[genome_1], fasta_dict[genome_2], fasta_dict[genome_1], fasta_dict[genome_2]]
        identical_bools = [False, False, True, True]

        for blast_file_name, contig_dict, bool in zip(files, contig_dicts, identical_bools):
            blast_rows_df = append_FASTA_neighborhood_contigs(blast_file_name, AMR_gene, contig_dict, identical=bool)
            df = df.append(blast_rows_df)
            df.reset_index(drop=True, inplace=True)

        AMR_dict[genome_1 + '_' + genome_2 + '.blast.txt'] = df

    return AMR_dict


def get_blast_dict(blast_folder_paths):
    """
    Creates dictionary for all AMR genes where each keys is an AMR gene name (str) and each value is the dictionary
    containing genome combinations that were blasted together as keys and their accompanying blast file data in a
    dataframe as values.
    """
    BLAST_df_dict = {}
    for blast_amr_subdir in blast_folder_paths:
        # For each gene, get dictionary of all blast results above %ID criterion
        blast_gene_dict = make_BLAST_df_dict(blast_amr_subdir)

        # Store data
        AMR_gene = get_filename(blast_amr_subdir).split('.')[0]
        BLAST_df_dict[AMR_gene] = blast_gene_dict

    return BLAST_df_dict


def make_BLAST_df_dict(blast_gene_dir):
    """
    Returns a dictionary of BLAST dataframes, where a dataframe exists for each BLAST file and the name is the key.
    """
    blast_dataframes = {}
    if os.path.isdir(blast_gene_dir):
        for filename in os.listdir(blast_gene_dir):
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
    Retains only protein sequences that share at least 70% overall sequence identity (i.e., homologs) in BLAST dict.
    """
    filtered_BLAST_dict = {}
    for blast_file, blast_df in blast_gene_dict.items():

        # Filter rows to drop all rows with percent identity lower than specified threshold
        df = blast_df[blast_df['PI'] >= cutoff]
        df.reset_index(drop=True, inplace=True)

        # Remove rows where query ID and sub ID are the same: was only needed for normalized bitscores calculation
        final_df = df[df['query_id'] != df['sub_id']]
        final_df.reset_index(drop=True, inplace=True)

        filtered_BLAST_dict[blast_file] = df

    return filtered_BLAST_dict


def remove_BLAST_duplicates(blast_df):
    """
    Given a BLAST file dataframe, removes duplicate entries for same gene by keeping only the highest bitscore.
    """
    blast_df = blast_df.sort_values('normalized_bitscore', ascending=False)
    blast_df = blast_df.drop_duplicates(['query_id'])
    blast_df.reset_index(drop=True, inplace=True)

    return blast_df


def get_neighborhood_from_fasta(fasta_file_path):
    """
    Given a FASTA file of an AMR gene neighborhood created using extraction.py, creates a list of all the gene
    identifiers in the neighborhood without their protein sequences.
    """
    neighborhood = []
    with open(fasta_file_path, 'r') as infile:
        for line in infile.readlines():
            if line.startswith('>'):
                gene_id = line.replace('>', '').replace('\n', '')
                neighborhood.append(gene_id)

    return neighborhood


def get_neighborhoods_dict(fasta_dir_path):
    """
    Given user-provided path to FASTA files for the genomes being analyzed, returns a dictionary that contains the
    neighborhood of that AMR gene in each genome for reference for later construction of similarity matrices.
    Keys: AMR genes
    Values: Genome dict
           - Keys: Genome ID (str)
           - Values: Neighborhood
    """
    neighborhoods = {}
    for amr_dir in get_full_filepaths(fasta_dir_path):
        AMR_gene = get_filename(amr_dir)
        amr_gene_neighborhoods = {}
        for fasta_filepath in get_full_filepaths(amr_dir):
            fasta_filename = get_filename(fasta_filepath)
            neighborhood = get_neighborhood_from_fasta(fasta_filepath)
            amr_gene_neighborhoods[fasta_filename] = neighborhood
        neighborhoods[AMR_gene] = amr_gene_neighborhoods

    return neighborhoods


def get_similarity_matrices(amr_blast_dict, neighborhoods_dict, neighborhood_size, fasta_path):
    """
    Calculates similarity matrix for each AMR gene using data from BLAST dataframes and extracted neighborhoods
    """
    similarity_matrices_dict = {}
    AMR_gene_genomes_dict = {}

    for AMR_gene, blast_df_dict in amr_blast_dict.items():

        # Get genome IDs for every genome that was used to generate the BLAST results
        genome_ids = [get_filename(filepath) for filepath in get_full_filepaths(fasta_path + '/' + AMR_gene)]

        # store in AMR gene genomes dict
        AMR_gene_genomes_dict[AMR_gene] = genome_ids

        # Set diagonal to max similarity score possible (for each neighborhood's score against itself)
        max_similarity_score = neighborhood_size * 2 + 1
        similarity_matrix = [[max_similarity_score for _ in range(len(genome_ids))] for _ in range(len(genome_ids))]

        for permutation in itertools.permutations(range(len(genome_ids)), 2):

            # Get genome indices and retrieve respective BLAST file data
            genome_1_index = permutation[0]
            genome_2_index = permutation[1]

            # Get genome names
            genome_1_id = genome_ids[genome_1_index]
            genome_2_id = genome_ids[genome_2_index]

            # Retrieve BLAST dataframe comparing the two genomes
            try:
                blast_filename = genome_1_id + '_' + genome_2_id + '.blast.txt'
                blast_df = blast_df_dict[blast_filename]

                # Get gene neighborhoods for comparison of neighborhood sizes (i.e. contig ends)
                neighborhood_1 = neighborhoods_dict[AMR_gene][genome_1_id]
                neighborhood_2 = neighborhoods_dict[AMR_gene][genome_2_id]

            except KeyError:
                try:
                    blast_filename = genome_2_id + '_' + genome_1_id + '.blast.txt'
                    blast_df = blast_df_dict[blast_filename]

                    # Get gene neighborhoods for comparison of neighborhood sizes (i.e. contig ends)
                    neighborhood_1 = neighborhoods_dict[AMR_gene][genome_2_id]
                    neighborhood_2 = neighborhoods_dict[AMR_gene][genome_1_id]
                    
                except KeyError:
                    print("Key error for gene {g}: {a} and {b} didn't exist...".format(g=AMR_gene,
                                                                                       a=genome_1_id + '_' + genome_2_id + '.blast.txt',
                                                                                       b=genome_2_id + '_' + genome_1_id + '.blast.txt'))

            # ORIGINAL FILTERING CONDITION
            blast_neighborhood = blast_df[(((blast_df['query_id'].isin(neighborhood_1)) &
                                            (blast_df['sub_id'].isin(neighborhood_2))) |
                                           (((blast_df['query_id'].isin(neighborhood_2)) &
                                             (blast_df['sub_id'].isin(neighborhood_1)))))]

            blast_neighborhood_df = remove_BLAST_duplicates(blast_neighborhood)

            sum_1 = round(blast_neighborhood_df['normalized_bitscore'].sum(), 3)

            if len(neighborhood_1) > len(neighborhood_2):
                len_diff = len(neighborhood_1) - len(blast_neighborhood_df)
                sum_2 = sum_1 + len(neighborhood_1) - len_diff - len(blast_neighborhood_df)

            elif len(neighborhood_2) > len(neighborhood_1):
                len_diff = len(neighborhood_2) - len(blast_neighborhood_df)
                sum_2 = sum_1 + len(neighborhood_2) - len_diff - len(blast_neighborhood_df)

            else:
                sum_2 = sum_1

            # Update values in upper and lower diagonal
            final_sum = round(sum_2, 3)
            similarity_matrix[genome_1_index][genome_2_index] = final_sum
            similarity_matrix[genome_2_index][genome_1_index] = final_sum

        similarity_matrices_dict[AMR_gene] = similarity_matrix

    return similarity_matrices_dict, AMR_gene_genomes_dict


def get_average_similarity_score(np_similarity_matrix, genome_names):
    """
    Computes average similarity score across all genomes for a given AMR gene's similarity matrix.
    """
    df = pd.DataFrame(data=np_similarity_matrix, index=genome_names, columns=genome_names)
    values = df.values

    return int(round(values.mean()))


def get_distance_matrix(np_similarity_matrix, genome_names):
    """"
    Converts numpy symmetric matrix array into a distance matrix dataframe!
    """
    df = pd.DataFrame(data=np_similarity_matrix, index=genome_names, columns=genome_names)
    values = df.values
    df = df.transform(lambda x: 1 - x / values.max())
    df = df.round(decimals=3)

    return df


def get_maximum_distance_score(distance_matrix, genome_names):
    df = pd.DataFrame(data=distance_matrix, index=genome_names, columns=genome_names)
    values = df.values

    return values.mean()


def get_sparse_matrix(np_matrix):
    """
    Converts a numpy distance matrix to a scipy sparse matrix (important for use with MCL clustering).
    """
    sparse_matrix = sparse.csr_matrix(np_matrix)

    return sparse_matrix

def get_newick_from_tree(tree, ):
    return

def UPGMA_clustering(condensed_distance_matrix):
    """
    Applies UPGMA clustering to neighborhood similarity or symmetric distance matrix.
    """
    try:
        upgma_linkage = hierarchy.average(condensed_distance_matrix)
        return upgma_linkage
    except ValueError:
        print("Empty distance matrix was passed for the gene!")


#def MCL_hyperparameter_tuning(sparse_distance_matrix):
#    """
#    Performs hyperparameter tuning to find optimal inflation parameter for a given sparse distance matrix we want to
#    apply MCL to using an Optuna study.
#    """
#    def objective(trial):
#        """
#        Optuna optimization trial for Markov Clustering Algorithm modularity (Q) score
#        (see documentation at: https://markov-clustering.readthedocs.io/en/latest/readme.html#choosing-hyperparameters).
#        """
#        inflation = trial.suggest_float('inflation', 1.0, 15.0)

#        result = mcl.run_mcl(sparse_distance_matrix, inflation=inflation)
#        clusters = mcl.get_clusters(result)

#        Q_score = mcl.modularity(matrix=result, clusters=clusters)
#        return Q_score

#    study = optuna.create_study(
#        study_name='Markov_Clustering_Optimization',
#        direction='maximize',
#        pruner=optuna.pruners.HyperbandPruner(max_resource='auto')
#    )
#    study.optimize(objective, n_trials=100)

    # Output results
 #   print("Best Q obtained during optimization: ", study.best_value)
 #   print("Best inflation parameter: ", study.best_params)

 #   return study.best_params.get('inflation')


def MCL_clustering(matrix, inflation):
    """
    Performs MCL clustering on a neighborhood sparse similarity or symmetric distance matrix.
    """
    sparse_matrix = get_sparse_matrix(matrix)
    #inflation = MCL_hyperparameter_tuning(sparse_matrix)
    result = mcl.run_mcl(sparse_matrix, inflation=inflation)
    clusters = mcl.get_clusters(result)

    return clusters


#def DBSCAN_hyperparameter_tuning(np_distance_matrix):
#    """
#    Performs hyperparameter tuning to find optimal minPts and epsilon values for a given numpy distance matrix
#    we want to apply DBSCAN to using an Optuna study.
#    """
#
#    def objective(trial):
#        """
#        Optuna optimization trial for DBSCAN Density Based Validation (DBCV) score
#        (see documentation at: https://markov-clustering.readthedocs.io/en/latest/readme.html#choosing-hyperparameters).
#        """
#        epsilon = trial.suggest_float('eps', 0, 1)
#        min_points = trial.suggest_int('min_samples', 1, 5)

#        distance_matrix = StandardScaler().fit_transform(np_distance_matrix)
#        dbscan = DBSCAN(eps=epsilon, min_samples=min_points).fit(distance_matrix)
#        dbcv_score = DBCV(X=distance_matrix, labels=dbscan.labels_, dist_function=euclidean)

#        return dbcv_score

#    study = optuna.create_study(
#        study_name='DBSCAN_Optimization',
#        direction='maximize',
#        pruner=optuna.pruners.HyperbandPruner(max_resource='auto')
#    )
#    study.optimize(objective, n_trials=100)

    # Output results
#    print("Best DBCV score obtained during optimization: ", study.best_value)
#    print("Best epsilon, min_samples parameter values found: ", study.best_params)

#    return study.best_params.get('eps'), study.best_params.get('min_samples')


def DBSCAN_clustering(np_distance_matrix, epsilon, minpts):
    """
    Applies DBSCAN clustering to a neighborhood similarity or symmetric distance matrix.
    """
    #epsilon, min_points = DBSCAN_hyperparameter_tuning(np_distance_matrix)
    distance_matrix = StandardScaler().fit_transform(np_distance_matrix)
    dbscan = DBSCAN(eps=epsilon, min_samples=minpts).fit(distance_matrix)

    return dbscan, dbscan.labels_


def check_clustering_savepaths(output_path):
    """
    Ensures that within the specified output directory, there is a directory for clustering visualizations with
    a distinct subdirectory containing the clustering results graph for each AMR gene with that algorithm.
    """
    for clustering_type in ['UPGMA', 'DBSCAN', 'MCL']:
        save_path = output_path + '/clustering/' + clustering_type
        check_output_path(save_path)


def load_JSON_data(output_path, AMR_gene, surrogates=False):
    """
    Helper function for loading JSON neighborhood data.
    """
    json_data = ''
    if surrogates:
        gene_path = output_path + '/JSON/' + AMR_gene + '_surrogates.json'
    else:
        gene_path = output_path + '/JSON/' + AMR_gene + '.json'

    print(gene_path)
    with open(gene_path, 'r') as infile:
        if len(infile.readlines()) != 0:
            infile.seek(0)
            json_data = json.load(infile)

    return json_data, gene_path


def update_JSON_links_PI(BLAST_df_dict, output_path, surrogates=False):
    """
    Updates JSON representations of AMR gene neighborhoods created using extraction module so that gene cluster links
    reflect percent identities found in blast results.
    """
    for AMR_gene, blast_files_dict in BLAST_df_dict.items():

        # Load AMR gene JSON link data
        json_data, gene_path = load_JSON_data(output_path, AMR_gene, surrogates)

        # Update each link according to the respective blast results
        for i in range(len(json_data["links"])):

            # Contig identifiers
            contig_data = json_data["links"][i]["uid"].split('-')
            contig_id_1 = contig_data[0].split('_')
            contig_id_2 = contig_data[1].split('_')
            contig_1 = contig_id_1[1] + '_' + contig_id_1[2]
            contig_2 = contig_id_2[1] + '_' + contig_id_2[2]

            # Genome names
            genome_1 = contig_1.split('_')[0]
            genome_2 = contig_2.split('_')[0]

            # Check respective BLAST file dataframe and update percent identity
            try:
                df = BLAST_df_dict[AMR_gene][genome_1 + '_' + genome_2 + '.blast.txt']
                row = df.loc[((df['query_id'] == contig_1) & (df['sub_id'] == contig_2))]
                PI = row.PI.tolist()[0]

            except KeyError:
                try:
                    df = BLAST_df_dict[AMR_gene][genome_2 + '_' + genome_1 + '.blast.txt']
                    row = df.loc[((df['query_id'] == contig_1) & (df['sub_id'] == contig_2))]
                    PI = row.PI.tolist()[0]
                except KeyError:
                    PI = 0.70
                except IndexError:
                    PI = 0.70

            except IndexError:
                PI = 0.70

            try:
                df = BLAST_df_dict[AMR_gene][genome_2 + '_' + genome_1 + '.blast.txt']
                row = df.loc[((df['query_id'] == contig_1) & (df['sub_id'] == contig_2))]
                PI = row.PI.tolist()[0]

            except KeyError:
                try:
                    df = BLAST_df_dict[AMR_gene][genome_1 + '_' + genome_2 + '.blast.txt']
                    row = df.loc[((df['query_id'] == contig_1) & (df['sub_id'] == contig_2))]
                    PI = row.PI.tolist()[0]
                except KeyError:
                    PI = 0.70
                except IndexError:
                    PI = 0.70

            except IndexError:
                PI = 0.70

            json_data["links"][i]["identity"] = PI

        # Overwrite JSON file with updated data
        with open(gene_path, 'w') as outfile:
            json.dump(json_data, outfile)


def map_genome_id_to_dendrogram_leaves(upgma_clusters, genome_to_num_mapping):
    """
    Given a UPGMA dendrogram and a dictionary mapping genomes to integers corresponding to the indices of the distance
    matrix used to generate the dendrogram, returns a list of genome names in order of the dendrogram leaves from left
    to right.
    """
    upgma_leaves = upgma_clusters['ivl']
    genome_order = []
    for leaf in upgma_leaves:
        genome = genome_to_num_mapping[leaf]
        genome_order.append(genome)

    return genome_order


def order_cluster_data_by_dendrogram(genome_order_dict, json_cluster_data):
    """
    Given JSON cluster data in the clustermap format, reorders cluster data according to dendrogram leaves from
    left to right.
    """
    clusters = []
    for genome in genome_order_dict:
        for cluster in json_cluster_data:
            if cluster["name"] == genome:
                print(cluster)
                clusters.append(cluster)
    return clusters


def order_JSON_clusters_UPGMA(output_path, AMR_gene, upgma_clusters, genome_to_num_mapping, surrogates=False):
    """
    Reorders how genomes are encoded in their respective JSON files for an AMR gene according to how they
    were clustered by UPGMA (from left to right).
    """
    # Load AMR gene JSON cluster data
    json_data = ''
    if surrogates:
        gene_path = output_path + '/JSON/' + AMR_gene + '_surrogates.json'
    else:
        gene_path = output_path + '/JSON/' + AMR_gene + '.json'

    with open(gene_path, 'r') as infile:
        if len(infile.readlines()) != 0:
            infile.seek(0)
            json_data = json.load(infile)

    # Reorder cluster data genomes according to UPGMA leaves ordering
    genome_order = map_genome_id_to_dendrogram_leaves(upgma_clusters, genome_to_num_mapping)
    clusters = order_cluster_data_by_dendrogram(genome_order, json_data["clusters"])

    # Update JSON data
    json_data["clusters"] = clusters

    # Update file
    with open(gene_path, 'w') as outfile:
        json.dump(json_data, outfile)


def make_representative_UPGMA_cluster_JSON(output_path, AMR_gene, upgma_clusters, genome_to_num_mapping):
    """
    Creates a JSON file with one representative genome from each UPGMA cluster.
    """
    # Load AMR gene JSON cluster data
    json_data = ''
    with open(output_path + '/JSON/' + AMR_gene + '.json', 'r') as infile:
        if len(infile.readlines()) != 0:
            infile.seek(0)
            json_data = json.load(infile)

    # Determine genomes to keep in representation
    representative_cluster_genomes = []

    genome_order = map_genome_id_to_dendrogram_leaves(upgma_clusters, genome_to_num_mapping)
    cluster_df = pd.DataFrame({'genome': genome_order, 'cluster': upgma_clusters['leaves_color_list']})
    print(AMR_gene)
    print(cluster_df)

    unique_clusters = set(upgma_clusters['leaves_color_list'])
    for cluster in unique_clusters:
        cluster_rows = cluster_df.loc[cluster_df['cluster'] == cluster]
        representative_cluster_genomes.append(cluster_rows.iloc[0].genome)
    del cluster_df

    # Update cluster data to only include those genomes
    clusters = []
    for genome in representative_cluster_genomes:
        for cluster in json_data["clusters"]:
            if cluster["name"] == genome:
                print(cluster)
                clusters.append(cluster)

    # Update JSON data
    json_data["clusters"] = clusters

    #upgma_json_data = remove_defunct_clustermap_data(json_data)

    #print("UPGMA JSON DATA:")
    #print(upgma_json_data)

    # Update file
    with open(output_path + '/JSON/' + AMR_gene + '_upgma.json', 'w') as outfile:
        json.dump(json_data, outfile)

    # Make respective HTML file for Coeus
    write_clustermap_JSON_HTML(AMR_gene, '../sample_data', output_path, rep_type='upgma')


def cluster_neighborhoods(assembly_path, fasta_path, blast_path, output_path,
                          neighborhood_size=10, inflation=2, epsilon=0.5, minpts=5):
    """
    Driver script for clustering neighborhoods obtained from extraction module.
    """
    # If BLAST files do not exist, make them prior to clustering
    check_output_path(blast_path)
    blast_dir_files = os.listdir(blast_path)
    fasta_dir_files = os.listdir(fasta_path)

    fa_dir_files = os.listdir(assembly_path)
    num_genome_combinations = list(itertools.permutations(fa_dir_files, 2))
    #if len(blast_dir_files) < (len(num_genome_combinations) + len(fa_dir_files)):
    #    print("Creating whole genome DIAMOND All-vs-All BLAST files...")
    #    blast_neighborhoods(assembly_path, blast_path)

    #if len(blast_dir_files) < len(fasta_dir_files):
    #    blast_neighborhoods(assembly_path, blast_path)

    # Make BLAST dataframe dictionary
    print("Fetching relevant BLAST data from DIAMOND outputs for each respective neighborhood...")
    #BLAST_df_dict = get_blast_dict_whole_genomes(fasta_path, blast_path)

    blast_folder_paths = get_full_filepaths(blast_path)

    BLAST_df_dict = {}
    for blast_amr_subdir in blast_folder_paths:

        # For each gene, get dictionary of all blast results above %ID criterion
        blast_gene_dict = make_BLAST_df_dict(blast_amr_subdir)

        # Store data
        AMR_gene = get_filename(blast_amr_subdir).split('.')[0]
        BLAST_df_dict[AMR_gene] = blast_gene_dict

    # Calculate normalized bitscores and save as column
    print("Calculating normalized bitscores for downstream scoring...")
    final_BLAST_dict = {}
    for amr_gene_subdir, df_dict in BLAST_df_dict.items():
        blast_files_dict = {}
        for filename, blast_df in df_dict.items():
            print("Getting normalized bitscores for gene {a} in genome {f}...".format(f=filename, a=amr_gene_subdir))
            final_df = get_normalized_bitscores(blast_df)
            blast_files_dict[filename] = final_df

        # Store data
        AMR_gene = get_filename(amr_gene_subdir).split('.')[0]
        final_BLAST_dict[AMR_gene] = blast_files_dict
        #final_BLAST_dict[AMR_gene] = filter_BLAST_results(blast_files_dict)

    # Update UID genes in JSON neighborhood representation
    #update_unidentified_genes_data(BLAST_df_dict, output_path, surrogates=False)

    # Update links in each AMR gene JSON file to reflect percent identities of blast hits
    print("Updating JSON neighborhood representations' links percent identities according to BLAST results...")
    update_JSON_links_PI(BLAST_df_dict, output_path, surrogates=False)
    update_JSON_links_PI(BLAST_df_dict, output_path, surrogates=True)

    # Update links in each AMR gene JSON file to reflect percent identities of blast hits
    update_JSON_links_PI(BLAST_df_dict, output_path)

    # Get neighborhoods dict for calculating similarity matrices (needed to compare contig ends)
    for AMR_gene_subdir in os.listdir(fasta_path):
        # Remove non-FASTA files
        remove_files(fasta_path + '/' + AMR_gene_subdir, '.fasta')

    neighborhoods = get_neighborhoods_dict(fasta_path)

    # Calculate similarity matrix for each gene
    print("Calculating similarity matrices for each AMR gene...")
    check_output_path(output_path)
    similarity_matrices_dict, AMR_genome_names_dict = get_similarity_matrices(final_BLAST_dict,
                                                                              neighborhoods,
                                                                              neighborhood_size,
                                                                              fasta_path)
    # Make output dir for clustering results
    check_output_path(output_path + '/clustering')

    # Generate distance matrices
    print("Calculating distance matrices...")

    distance_matrices_df_dict = {}
    average_similarity_scores_dict = {}
    for AMR_gene, similarity_matrix_data in similarity_matrices_dict.items():

        similarity_matrix = np.array(similarity_matrix_data)
        genome_names = AMR_genome_names_dict[AMR_gene]

        try:
            distance_matrix = get_distance_matrix(similarity_matrix, genome_names)
            distance_matrices_df_dict[AMR_gene] = distance_matrix
            average_similarity_scores_dict[AMR_gene] = get_average_similarity_score(similarity_matrix, genome_names)

        # if directory empty of BLAST results
        except ValueError:
            print("Was unable to locate BLAST results for gene {g}. Skipped.".format(g=AMR_gene))
            pass

    print("Clustering neighborhoods...")
    check_clustering_savepaths(output_path)

    max_distance_scores_dict = {}

    for AMR_gene, distance_matrix_df in distance_matrices_df_dict.items():

        # Get distance matrix
        np.fill_diagonal(distance_matrix_df.values, 0)
        distance_matrix = np.array(distance_matrix_df.values)

        genome_names = AMR_genome_names_dict[AMR_gene]

        # Retain maximum distance score for each AMR gene for distances histogram
        max_distance_scores_dict[AMR_gene] = get_maximum_distance_score(distance_matrix, genome_names)

        # Reorganize JSON clusters according to UPGMA order (left to right)
        genome_to_num_mapping = {}
        for i in range(len(genome_names)):
            genome_to_num_mapping[str(i)] = genome_names[i]

        print("Generating UPGMA clusters for {g}...".format(g=AMR_gene))
        try:
            # Get UPGMA linkage matrix
            upgma_linkage = UPGMA_clustering(distance_matrix)

            # Use linkage matrix to build dendrogram used for updating JSON data: reorder clusters, make UPGMA view
            upgma_dendrogram = hierarchy.dendrogram(upgma_linkage)
            order_JSON_clusters_UPGMA(output_path, AMR_gene, upgma_dendrogram, genome_to_num_mapping, surrogates=False)
            order_JSON_clusters_UPGMA(output_path, AMR_gene, upgma_dendrogram, genome_to_num_mapping, surrogates=True)
            make_representative_UPGMA_cluster_JSON(output_path, AMR_gene, upgma_dendrogram, genome_to_num_mapping)

            # Interactive dendrogram visualization
            condensed_distance_matrix = squareform(distance_matrix)
            print(condensed_distance_matrix)
            #plotly_dendrogram(condensed_distance_matrix, genome_names, AMR_gene, output_path)
            #plotly_dendrogram(upgma_linkage, genome_names, AMR_gene, output_path)
            plotly_dendrogram(distance_matrix, genome_names, AMR_gene, output_path)

            # Save distance matrix as textfile
            check_output_path(output_path + '/clustering/distance_matrices')
            distance_matrix_df.to_csv(output_path + '/clustering/distance_matrices/' + AMR_gene + '_distance_matrix.csv', sep='\t', index=False)

        except IndexError:
            print("Unable to perform UPGMA clustering for gene {g}. " \
                  "UPGMA results for {g} will be omitted.".format(g=AMR_gene))

        except TypeError:
            print("Unable to perform UPMA clustering for gene {g}. " \
                  "UPGMA results for {g} will be omitted.".format(g=AMR_gene))

        print("Generating DBSCAN clusters for {g}...".format(g=AMR_gene))

        try:
            dbscan_clusters, labels = DBSCAN_clustering(distance_matrix, epsilon, minpts)

            # Plot DBSCAN clusters using distance matrix and PCoA: colour according to cluster assignment
            plotly_pcoa(distance_matrix_df, genome_names, labels, AMR_gene, output_path)

        except IndexError:
            print("Unable to perform DBSCAN clustering for gene {g}. " \
                  "DBSCAN results for {g} will be omitted.".format(g=AMR_gene))

        except KeyError:
            print("Unable to perform DBSCAN clustering for gene {g}. " \
                  "DBSCAN results for {g} will be omitted.".format(g=AMR_gene))

    for AMR_gene, similarity_matrix_data in similarity_matrices_dict.items():

        try:
            similarity_matrix = np.array(similarity_matrix_data)
            genome_names = AMR_genome_names_dict[AMR_gene]

            df = pd.DataFrame(data=similarity_matrix, index=genome_names, columns=genome_names)
            values = df.values

            print("Generating Markov clusters for {g}...".format(g=AMR_gene))

            clusters = MCL_clustering(similarity_matrix, inflation)
            sparse_sim_matrix = get_sparse_matrix(similarity_matrix)
            plotly_mcl_network(sparse_sim_matrix, clusters, genome_names, AMR_gene, output_path)

            # Save distance matrix as textfile
            check_output_path(output_path + '/clustering/similarity_matrices')
            similarity_matrix_df = pd.DataFrame(data=similarity_matrix, index=genome_names, columns=genome_names)
            similarity_matrix_df.to_csv(output_path + '/clustering/similarity_matrices/' + AMR_gene + \
                                          '_similarity_matrix.csv', sep='\t', index=False)

        except IndexError:
            print("Unable to perform MCL clustering for gene {g}. "\
                  "MCL results for {g} will be omitted.".format(g=AMR_gene))

        except ValueError:
            print("Unable to perform MCL clustering for gene {g}. "\
                  "MCL results for {g} will be omitted.".format(g=AMR_gene))

    # Generate summary histograms for analyzed genomes' similarities
    print("Generating average similarity and max distance histograms for all neighborhoods...")
    plot_similarity_histogram(average_similarity_scores_dict, output_path)
    plot_distance_histogram(max_distance_scores_dict, output_path)


def main(args=None):
    args = parse_args(args)
    cluster_neighborhoods(args.ASSEMBLY_PATH, args.FASTA_PATH, args.BLAST_PATH, args.OUTPUT_PATH, args.n, args.i, args.e, args.m)


if __name__ == '__main__':
    sys.exit(main())
