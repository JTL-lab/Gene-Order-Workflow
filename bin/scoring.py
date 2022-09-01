#!/usr/bin/env python

"""
Script for scoring utilities.
"""

def calculate_bitscore(g1_g1_bitscore, g1_g2_bitscore):
    return round(float(g1_g1_bitscore / g1_g2_bitscore), 3)

def get_normalized_bitscores(blast_df):
    """
    Adds normalized bitscore column to BLAST results dataframe given dataframe contains data from BLAST textfiles
    in output format 6.
    """
    blast_df['normalized_bitscore'] = blast_df.bitscore
    blast_df_copy = blast_df.copy()

    # Get dataframe containing only rows where the query and comparison sequences are identical
    #identical_df = blast_df[blast_df['query_id'] == blast_df['sub_id']]
    identical_df = blast_df_copy.query('query_id == sub_id')
    identical_df.reset_index(drop=True, inplace=True)

    # Get dataframe containing only rows where query and comparison sequences are different
    #unique_df = blast_df[blast_df['query_id'] != blast_df['sub_id']]
    unique_df = blast_df.query('query_id != sub_id')
    unique_df.reset_index(drop=True, inplace=True)

    for index, gene_1 in identical_df.iterrows():

        # Get the corresponding row in the BLAST dataframe and set its normalized bitscore to 1
        query = gene_1['query_id']
        blast_df.loc[(blast_df['query_id'] == query) & (blast_df['sub_id'] == query), 'normalized_bitscore'] = 1

        # Retain bitscore for ahead calculations
        g1_g1_bitscore = gene_1['bitscore']

        # Update normalized bitscore for each row without
        for index, gene_2 in unique_df.iterrows():
            if gene_2['query_id'] == query:

                # Obtain normalized bitscore
                g1_g2_bitscore = gene_2['bitscore']
                normalized_bitscore = calculate_bitscore(g1_g1_bitscore, g1_g2_bitscore)

                # Locate corresponding row in BLAST_df and update its normalized bitscore col value
                sub_query = gene_2['sub_id']
                blast_df.loc[(blast_df['query_id'] == query) & (blast_df['sub_id'] == sub_query), 'normalized_bitscore'] = normalized_bitscore

    return blast_df
