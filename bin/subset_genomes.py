#!/usr/bin/env python

"""
Script to subset genomes at pangenomic level. Given a similarity threshold P based on gene presence-absence from
pangenome analysis with Panaroo, retains only one copy of each pair of genomes sharing P% genes for downstream analysis.
"""

import argparse
import itertools
import csv
import sys
import os
import pandas as pd
from itertools import combinations
from utils import check_output_path
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import AgglomerativeClustering

def parse_args(args=None):
    Description = "Subsets genomes based on percent similarity of gene presence-absence profile from pangenomic analysis."
    Epilog = "Example usage: python subset_genomes.py <PANAROO_OUTPUT_PATH> <OUTPUT_PATH> -p <PERCENTAGE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('PANAROO_OUTPUT_PATH', metavar='pan_path', type=str,
                        help='Path to directory containing Panaroo outputs')
    parser.add_argument('OUTPUT_PATH', metavar='output_path', type=str, help='Path to output directory where '
                                                                             'CSV file will be outputted.')
    parser.add_argument('-p', metavar='p', type=float, default=0.99, help='Similarity threshold for genome comparison.')
    return parser.parse_args(args)


def get_binary_csv(panaroo_output_path):
    """Loads Panaroo presence-absence CSV and converts it to a presence-absence matrix indexed by gene names."""
    pa_df = pd.read_csv(panaroo_output_path + 'gene_presence_absence.csv', index_col=0)
    pa_matrix = pa_df.iloc[:, 3:].notnull().astype(int)
    return pa_matrix


def get_gene_counts(pa_matrix):
    """Given presence-absence matrix, returns a list of tuples for each gene where the first value is the gene name and
    the second value is the number of genomes it occurs in."""
    gene_counts = []
    for gene in pa_matrix.index:
        gene_counts.append((gene, pa_matrix.loc[gene].sum()))
    return gene_counts


def get_core_genes(gene_counts, num_genomes, core_threshold=0.99):
    """Returns a list of all core genes"""
    threshold = num_genomes * core_threshold
    core_genes = [gene for gene, count in gene_counts if count >= threshold]
    return core_genes


def group_genomes_by_pa(core_genes, pa_matrix, similarity_threshold):
    """Returns groupings of genomes based on shared percentage p% of core genome."""
    core_pa_matrix = pa_matrix.loc[core_genes]

    similarity_matrix = 1 - squareform(pdist(core_pa_matrix.T, metric='jaccard'))

    clustering = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='average',
                                         distance_threshold=1 - similarity_threshold)
    clustering.fit(1 - similarity_matrix)

    # Group genomes by clusters
    clusters = {i: [] for i in range(clustering.n_clusters_)}
    for genome, cluster_id in zip(core_pa_matrix.columns, clustering.labels_):
        clusters[cluster_id].append(genome)

    genome_clusters_dict = {}
    for cluster in clusters.values():
        if cluster:
            # Set the first genome as the key and the rest as the value
            if len(cluster) > 1:
                genome_clusters_dict[cluster[0]] = cluster[1:]
            else:
                genome_clusters_dict[cluster[0]] = []

    return genome_clusters_dict


def get_representative_genomes(genome_clusters):
    """Returns dict of representative genomes as keys and list of corresponding surrogate genomes as values"""
    representative_genomes = dict()
    for i in range(len(genome_clusters)):
        rep_genome = genome_clusters[i][0]
        representative_genomes[rep_genome] = genome_clusters[i].pop(rep_genome)
    return representative_genomes


def subset_genomes(panaroo_output_path, output_path, similarity_threshold):
    """Driver script to obtain list of distinct, differentiated genomes to retain for downstream analysis and their surrogate genomes."""

    pa_matrix = get_binary_csv(panaroo_output_path)
    gene_counts = get_gene_counts(pa_matrix)

    num_genomes = len(pa_matrix.columns) - 1
    core_genes = get_core_genes(gene_counts, num_genomes)

    representative_genomes = group_genomes_by_pa(core_genes, pa_matrix, similarity_threshold)
    with open(output_path.strip('/') + '/subsetted_genomes.txt', 'w') as outfile:
        for representative, surrogates in representative_genomes.items():
            outfile.write(f"{representative}: {surrogates}\n")


def main(args=None):
    args = parse_args(args)
    subset_genomes(args.PANAROO_OUTPUT_PATH, args.OUTPUT_PATH, args.p)


if __name__ == '__main__':
    sys.exit(main())
