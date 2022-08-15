#!/usr/bin/env python

"""
Utils for file manipulation and general data preprocessing used by multiple modules.
"""

import os
import glob

def strip_brackets(var):
    """
    Preprocesses list or string with brackets to be bracketless string (e.g. for easier manipulation of dataframes).
    """
    if isinstance(var, list):
        final = str(var)
        clean_str = final.strip('[]')
    else:
        clean_str = var.strip('[]')

    return clean_str

def get_full_filepaths(dir_path):

    filenames = os.listdir(dir_path)
    filepaths = [dir_path + '/' + filenames[i] for i in range(len(filenames))]

    return filepaths

def get_dir_filepaths(filepaths, file_ext):
    """
    Given a path to the folder where all genome AMR neighborhoods for a given AMR gene are stored, returns a list
    of all neighborhood fasta filepaths.
    """
    try:
        neighborhood_filepaths = glob.glob(os.path.join(filepaths, file_ext))
        return neighborhood_filepaths
    except FileNotFoundError:
        print("Error: no {} files were found in the specified directory.".format(file_ext))

def get_filename(filepath):
    """
    Gets filename from filepath (i.e. removes directories and file suffix).
    """
    filename = os.path.basename(filepath).split(".")[0]
    return filename

def get_dir_filenames(filepaths):
    """
    Given list of filepaths, returns their names (i.e. filename stripped of suffix).
    """
    filenames = []
    for file_path in filepaths:
        filenames.append(os.path.basename(file_path).split(".")[0])

    return filenames

def make_output_file(output_path):
    """
    Creates an empty textfile in the specified output directory.
    """
    with open(output_path, 'w') as file:
        pass
    file.close()

def move_file(new_path, current_path):
    """
    Moves file from src to dest directory.
    """
    full_new_path = os.path.join(current_path, new_path)
    os.rename(current_path, full_new_path)

def remove_files(path_to_dir, file_ext):
    """
    Removes all files from a specified directory that do not have the given extension.
    """
    for file in os.listdir(path_to_dir):
        filepath = os.path.join(path_to_dir, file)
        if not filepath.endswith(file_ext):
            os.remove(filepath)