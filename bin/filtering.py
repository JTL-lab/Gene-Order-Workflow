#!/usr/bin/env python

"""
Given RGI files and their corresponding GBK files for complete bacterial genomes, as well as gene neighborhood clusters
generated using neighborhood_clustering.py, the utilities in this class are used for two main types of filtering:

a) Filtering of redundant AMR gene neighborhoods based on the clustering results for cleaner visualization downstream.
b) If user opts to investigate Loose hits using predictive module, filtering of low percent identity Loose hits.
"""

def filter_neighborhoods():
    return

def filter_Loose_hits(RGI_dataframes_dict, PI=0.80):
    """
    Returns RGI dataframes with only Loose hits above or equal to a given threshold of Percent Identity.
    If user does not specify PI for Best_Identities as cutoff, the default value is 80%.
    """
    for genome_id, RGI_df in RGI_dataframes_dict.items():
        # Retrieve all Loose hits for the given genome
        try:
            filtered_df = RGI_df[~((RGI_df['Cut_Off'] == 'Loose') and (RGI_df['Best_Identities'] < PI))]
            RGI_dataframes_dict[genome_id] = filtered_df
        except KeyError:
            print("No Loose hits found in the genome!")

    return RGI_dataframes_dict

