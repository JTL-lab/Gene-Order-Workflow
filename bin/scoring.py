#!/usr/bin/env python

"""
Script for scoring utilities.
"""
import pandas as pd

def get_normalized_bitscores(blast_df):
    """
    Adds normalized bitscore column to BLAST results dataframe given dataframe contains data from BLAST textfiles
    in output format 6.
    """
    # Get dataframe copy to modify
    blast_df_copy = blast_df.copy()

    # Use dictionary to store g1_g1 bitscores for each respective unique query id
    blast_df_copy['bitscore'] = blast_df_copy['bitscore'].where(blast_df_copy['query_id']==blast_df_copy['sub_id'], other=2.000)
    g1_g1_bitscore_dict = dict()
    for query_id, bitscore in zip(blast_df_copy.query_id.to_list(), blast_df_copy.bitscore.to_list()):
        if bitscore != 2.000:
            g1_g1_bitscore_dict[query_id] = bitscore

    # Calculate normalized bitscore by performing g1_g2_bitscore / g1_g1_bitscore as row-wise operator
    blast_df['normalized_bitscore'] = blast_df_copy.apply(lambda row: row.bitscore / g1_g1_bitscore_dict[row.query_id], axis=1)

    return blast_df