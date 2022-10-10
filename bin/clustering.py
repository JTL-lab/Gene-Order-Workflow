#!/usr/bin/env python

"""
Functions to cluster neighborhoods identified using extraction.py.
"""

import argparse
import itertools
import os
import shutil
import optuna
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster import hierarchy
from scipy.spatial.distance import euclidean
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

import markov_clustering as mcl

from DBCV import DBCV
from scoring import get_normalized_bitscores
from utils import get_filename, check_output_path, get_full_filepaths, remove_files
from visualization import plot_similarity_histogram, plot_distance_histogram, \
                          graph_UPGMA_clusters, draw_MCL_graph, graph_DBSCAN_clusters, \
                          plotly_pcoa, plotly_dendrogram, plotly_mcl_network


def parse_args(args=None):
    Description = "Cluster extracted AMR gene neighborhoods to compare conservation characteristics across genomes."
    Epilog = "Example usage: python clustering.py <FASTA_PATH> <OUTPUT_PATH> -n <NEIGHBORHOOD_SIZE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
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

    return parser.parse_args(args)


def run_BLAST(neighborhoods_filepaths, root_path, dir_output_path):
    """
    Given filepaths for all genome neighborhoods of an AMR gene X, runs All-vs-All BLAST
    and stores the results in output/BLAST/X.
    """
    # BLAST each neighborhood against every other neighborhood once
    for neighborhood_fasta_1, neighborhood_fasta_2 in itertools.combinations(neighborhoods_filepaths, 2):
        print("Blasting All-vs-All for {a} vs {b}...".format(a=neighborhood_fasta_1, b=neighborhood_fasta_2))

        # Get the filepaths for the neighborhoods to BLAST this iteration
        fasta_1_filename = get_filename(neighborhood_fasta_1)
        fasta_2_filename = get_filename(neighborhood_fasta_2)

        # Blast each file against itself: need bitscore data for bitscore normalization later
        os.system('blastp -query ' + neighborhood_fasta_1 + ' -subject ' + neighborhood_fasta_1 + '-task blastp ' +
                  '-outfmt 6 -max_hsps 1 -out ' + fasta_1_filename + '.blast.txt')
        os.system('blastp -query ' + neighborhood_fasta_2 + ' -subject ' + neighborhood_fasta_2 + ' -task blastp ' +
                  '-outfmt 6 -max_hsps 1 -out ' + fasta_2_filename + '.blast.txt')

        blast_output_filename = fasta_1_filename + '_' + fasta_2_filename + ".blast.txt"

        # Make output dir path and switch directories
        blast_output_path = dir_output_path + '/' + blast_output_filename

        # Make BLAST db for the first neighborhood
        os.system('makeblastdb -in ' + neighborhood_fasta_1 + ' -parse_seqids -dbtype prot')

        # Run BLAST All-vs-All with specified e-value cutoff
        os.system('blastp -query ' + neighborhood_fasta_2 + ' -db ' + neighborhood_fasta_1 + \
                  ' -task blastp -outfmt 6 -max_hsps 1 -out ' + blast_output_filename)

        # Combine all three files into the BLAST All-vs-All of the two genomes against each other and delete originals
        read_files = [fasta_1_filename + '.blast.txt', fasta_2_filename + '.blast.txt']
        with open(blast_output_filename, 'a') as outfile:
            for file in read_files:
                with open(file, 'r') as infile:
                    outfile.write(infile.read())
                # Delete file after appending
                os.remove(file)

        # Move file to correct directory in at output_path
        curr_dir = os.getcwd()
        curr_blast_path = os.path.join(curr_dir, blast_output_filename)
        print("Moving {b} to {o}...".format(b=curr_blast_path, o=blast_output_path))
        shutil.move(curr_blast_path, blast_output_path)


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


def make_BLAST_genome_dict(AMR_dict, blast_df_dict):
    BLAST_genome_dict = {}
    for AMR_gene, data in AMR_dict.items():
        blast_df = blast_df_dict[AMR_gene]
        for blast_filename, blast_df in data.items():
            temp = []
            blast_df.reset_index(drop=True, inplace=True)

            for index in range(len(blast_df)):
                if blast_df['query_id'][index] in blast_df:
                    temp.append((blast_df['query_id'][index], blast_df['sub_id'][index],
                                 blast_df['PI'][index], blast_df['bitscore'][index]))
    return BLAST_genome_dict


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


def get_genome_id_from_blastfile(blast_filename):
    blast_tokens = blast_filename.split('_').split('.')
    return blast_tokens[0], blast_tokens[1]


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
                blast_filename = genome_2_id + '_' + genome_1_id + '.blast.txt'
                blast_df = blast_df_dict[blast_filename]

                # Get gene neighborhoods for comparison of neighborhood sizes (i.e. contig ends)
                neighborhood_1 = neighborhoods_dict[AMR_gene][genome_2_id]
                neighborhood_2 = neighborhoods_dict[AMR_gene][genome_1_id]

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


def sim_matrix_to_symmetric_matrix(similarity_matrix):
    """
    Uses the averages of scores from similarity matrix to convert it into a symmetric matrix for clustering.
    """
    symmetric_matrix = []
    for row in range(len(similarity_matrix)):
        average_score = 0
        for col in range(len(similarity_matrix)):
            if similarity_matrix[col][1] == similarity_matrix[row][0] and \
                    similarity_matrix[col][0] == similarity_matrix[row][1]:
                # Get rounded score approximated to 3 decimal points
                average_score = round((similarity_matrix[row][2] + similarity_matrix[col][2]) / 2, 3)

        symmetric_matrix.append((similarity_matrix[row][0], similarity_matrix[row][1], average_score))

    return symmetric_matrix


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


def UPGMA_clustering(np_distance_matrix):
    """
    Applies UPGMA clustering to neighborhood similarity or symmetric distance matrix.
    """
    # UPGMA clustering to obtain dendogram
    try:
        upgma_dendogram = hierarchy.average(np_distance_matrix)
        return upgma_dendogram
    except ValueError:
        print("Empty distance matrix was passed for the gene!")


def MCL_hyperparameter_tuning(sparse_distance_matrix):
    """
    Performs hyperparameter tuning to find optimal inflation parameter for a given sparse distance matrix we want to
    apply MCL to using an Optuna study.
    """

    def objective(trial):
        """
        Optuna optimization trial for Markov Clustering Algorithm modularity (Q) score
        (see documentation at: https://markov-clustering.readthedocs.io/en/latest/readme.html#choosing-hyperparameters).
        """
        inflation = trial.suggest_float('inflation', 1.0, 15.0)

        result = mcl.run_mcl(sparse_distance_matrix, inflation=inflation)
        clusters = mcl.get_clusters(result)
        Q_score = mcl.modularity(matrix=result, clusters=clusters)

        return Q_score

    study = optuna.create_study(
        study_name='Markov_Clustering_Optimization',
        direction='maximize',
        pruner=optuna.pruners.HyperbandPruner(max_resource='auto')
    )
    study.optimize(objective, n_trials=100)

    # Output results
    print("Best Q obtained during optimization: ", study.best_value)
    print("Best inflation parameter: ", study.best_params)

    return study.best_params.get('inflation')


def MCL_clustering(sparse_distance_matrix):
    """
    Performs MCL clustering on a neighborhood sparse similarity or symmetric distance matrix.
    """
    inflation = MCL_hyperparameter_tuning(sparse_distance_matrix)
    result = mcl.run_mcl(sparse_distance_matrix, inflation=inflation)
    clusters = mcl.get_clusters(result)

    return clusters


def DBSCAN_hyperparameter_tuning(np_distance_matrix):
    """
    Performs hyperparameter tuning to find optimal minPts and epsilon values for a given numpy distance matrix
    we want to apply DBSCAN to using an Optuna study.
    """

    def objective(trial):
        """
        Optuna optimization trial for DBSCAN Density Based Validation (DBCV) score
        (see documentation at: https://markov-clustering.readthedocs.io/en/latest/readme.html#choosing-hyperparameters).
        """
        epsilon = trial.suggest_float('eps', 0, 1)
        min_points = trial.suggest_int('min_samples', 1, 5)

        distance_matrix = StandardScaler().fit_transform(np_distance_matrix)
        dbscan = DBSCAN(eps=epsilon, min_samples=min_points).fit(distance_matrix)
        dbcv_score = DBCV(X=distance_matrix, labels=dbscan.labels_, dist_function=euclidean)

        return dbcv_score

    study = optuna.create_study(
        study_name='DBSCAN_Optimization',
        direction='maximize',
        pruner=optuna.pruners.HyperbandPruner(max_resource='auto')
    )
    study.optimize(objective, n_trials=100)

    # Output results
    print("Best DBCV score obtained during optimization: ", study.best_value)
    print("Best epsilon, min_samples parameter values found: ", study.best_params)

    return study.best_params.get('eps'), study.best_params.get('min_samples')


def DBSCAN_clustering(np_distance_matrix):
    """
    Applies DBSCAN clustering to a neighborhood similarity or symmetric distance matrix.
    """
    #epsilon, min_points = DBSCAN_hyperparameter_tuning(np_distance_matrix)
    distance_matrix = StandardScaler().fit_transform(np_distance_matrix)
    dbscan = DBSCAN(eps=0.9, min_samples=2).fit(distance_matrix)

    return dbscan, dbscan.labels_


def check_clustering_savepaths(output_path):
    """
    Ensures that within the specified output directory, there is a directory for clustering visualizations with
    a distinct subdirectory containing the clustering results graph for each AMR gene with that algorithm.
    """
    for clustering_type in ['UPGMA', 'DBSCAN', 'MCL']:
        save_path = output_path + '/clustering/' + clustering_type
        check_output_path(save_path)


def cluster_neighborhoods(fasta_path, blast_path, output_path, neighborhood_size=10):
    """
    Driver script for clustering neighborhoods obtained from extraction module.
    """
    # If BLAST files do not exist, make them prior to clustering
    check_output_path(blast_path)
    blast_dir_files = os.listdir(blast_path)
    if len(blast_dir_files) == 0:

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
                run_BLAST(full_genome_paths, root_path, dir_output_path)

    # Make BLAST dataframe dictionary
    blast_folder_paths = get_full_filepaths(blast_path)

    BLAST_df_dict = {}
    for blast_amr_subdir in blast_folder_paths:
        # For each gene, get dictionary of all blast results above %ID criterion
        blast_gene_dict = make_BLAST_df_dict(blast_amr_subdir)

        # Store data
        AMR_gene = get_filename(blast_amr_subdir).split('.')[0]
        BLAST_df_dict[AMR_gene] = blast_gene_dict

    # Calculate normalized bitscores and save as column
    final_BLAST_dict = {}
    for amr_gene_subdir, df_dict in BLAST_df_dict.items():
        blast_files_dict = {}
        for filename, blast_df in df_dict.items():
            print("Getting normalized bitscores for gene {a} in genome {f}...".format(f=filename, a=amr_gene_subdir))
            final_df = get_normalized_bitscores(blast_df)
            blast_files_dict[filename] = final_df

        # Store data
        AMR_gene = get_filename(amr_gene_subdir).split('.')[0]
        final_BLAST_dict[AMR_gene] = filter_BLAST_results(blast_files_dict)

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
            pass

    print("Performing clustering...")
    check_clustering_savepaths(output_path)

    max_distance_scores_dict = {}

    for AMR_gene, distance_matrix_df in distance_matrices_df_dict.items():

        # Get distance matrix
        np.fill_diagonal(distance_matrix_df.values, 0)
        distance_matrix = np.array(distance_matrix_df.values)

        genome_names = AMR_genome_names_dict[AMR_gene]

        # Retain maximum distance score for each AMR gene for distances histogram
        max_distance_scores_dict[AMR_gene] = get_maximum_distance_score(distance_matrix, genome_names)

        print("Generating UPGMA clusters for {g}...".format(g=AMR_gene))
        upgma_clusters = UPGMA_clustering(distance_matrix)
        graph_UPGMA_clusters(upgma_clusters, genome_names, AMR_gene, output_path)
        plotly_dendrogram(distance_matrix, genome_names, AMR_gene, output_path)

        print("Generating DBSCAN clusters...")
        dbscan_clusters, labels = DBSCAN_clustering(distance_matrix)

        # Plot DBSCAN clusters using distance matrix and PCoA: colour according to cluster assignment
        graph_DBSCAN_clusters(distance_matrix_df, dbscan_clusters, labels, AMR_gene, output_path)
        plotly_pcoa(distance_matrix_df, genome_names, labels, AMR_gene, output_path)

    for AMR_gene, similarity_matrix_data in similarity_matrices_dict.items():

        similarity_matrix = np.array(similarity_matrix_data)
        genome_names = AMR_genome_names_dict[AMR_gene]

        df = pd.DataFrame(data=similarity_matrix, index=genome_names, columns=genome_names)
        values = df.values

        print("Generating Markov clusters for {g}...".format(g=AMR_gene))
        sparse_sim_matrix = get_sparse_matrix(similarity_matrix)
        clusters = MCL_clustering(sparse_sim_matrix)
        draw_MCL_graph(sparse_sim_matrix, clusters,
                       edge_color="silver",
                       save_name=AMR_gene,
                       save_path=output_path)
        plotly_mcl_network(sparse_sim_matrix, clusters, genome_names, AMR_gene, output_path)

    # Generate summary histograms for analyzed genomes' similarities
    print("Generating average similarity and max distance histograms for all neighborhoods...")
    plot_similarity_histogram(average_similarity_scores_dict, output_path)
    plot_distance_histogram(max_distance_scores_dict, output_path)


def main(args=None):
    args = parse_args(args)
    cluster_neighborhoods(args.FASTA_PATH, args.BLAST_PATH, args.OUTPUT_PATH, args.n)


if __name__ == '__main__':
    sys.exit(main())
