#!/usr/bin/env python

"""
Script for scoring utilities.
"""

def get_normalized_bitscores(blast_df):

    # Remove all rows where query and comparison sequence are identical
    df = blast_df[blast_df['query_id'] == blast_df['sub_id']]
    df.reset_index(drop=True, inplace=True)

    query_dict = {}
    for row in range(len(df)):
        query_dict[df['query_id'][row]] = df['bitscore'][row]

    normalized_bitscores = []
    for row in range(len(blast_df)):
        normalized_bitscores.append(round(float(blast_df['bitscore'][row]) / query_dict[df['query_id'][row]], 3))

    blast_df['normalized_bitscore'] = normalized_bitscores

    return blast_df

