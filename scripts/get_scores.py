#!/usr/bin/env python

"""
Script for scoring utilities.
"""

def normalize_BLAST_bitscore(BLAST_bitscore, BLAST_query_ID):
    return round((float(BLAST_bitscore) / BLAST_query_ID), 3)

