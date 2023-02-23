#!/usr/bin/env python

"""
Script for scoring utilities.
"""
import pandas as pd
import cProfile
import io
import pstats

def profile(fnc):
    """
    Decorator that uses cProfile to profile a function
    """
    def inner(*args, **kwargs):
        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval

    return inner


def calculate_bitscore(g1_g1_bitscore, g1_g2_row):
    return float(round(g1_g2_row.bitscore / g1_g1_bitscore, 3))


#@profile
def get_normalized_bitscores(blast_df):
    """
    Adds normalized bitscore column to BLAST results dataframe given dataframe contains data from BLAST textfiles
    in output format 6.
    """
    # Get dataframe containing only rows where the query and comparison sequences are identical
    blast_df_copy = blast_df.copy()

    # Get dataframe containing only rows where query and comparison sequences are different
    norm_bitscores = []
    for row in blast_df.itertuples():
        g1_g1_row = blast_df_copy[((blast_df_copy['query_id'] == row.query_id) & (blast_df_copy['sub_id'] == row.query_id))]
        g1_g1_bitscore_row = g1_g1_row['bitscore'].tolist()
        g1_g1_bitscore = g1_g1_bitscore_row[0]
        norm_bitscore = calculate_bitscore(g1_g1_bitscore, row)
        norm_bitscores.append(norm_bitscore)

    blast_df['normalized_bitscore'] = norm_bitscores

    return blast_df
