#!/usr/bin/env python

"""
Utils for file manipulation and general data preprocessing used by multiple modules.
"""

import os
import glob
import string
import random

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
        
def get_fasta_filename(filepath):
    """
    Gets filename from filepath (i.e. removes directories and file suffix).
    """
    filename = os.path.basename(filepath).split(".fasta")[0]
    return filename


def get_filename(filepath, rgi=False):
    """
    Gets filename from filepath (i.e. removes directories and file suffix).
    """
    if rgi:
        if '.fna' in filepath:
            filename = os.path.basename(filepath).split('.')[0]
        else:
            filename = os.path.basename(filepath).split('_')[0]
    else:
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


def check_output_path(path):
    """
    Checks for presence of specified output path and creates the directory if it does not exist
    """
    if not os.path.exists(path):
        os.makedirs(path)


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
        if not file.endswith(file_ext):
            os.remove(path_to_dir + '/' + file)


def generate_alphanumeric_string(length):
    """
    Generate a random alphanumeric string of fixed length
    """
    str = string.ascii_lowercase
    return ''.join(random.choice(str) for i in range(length))


def write_clustermap_JSON_HTML(AMR_gene, sample_data_path, out_path, rep_type='standard'):

    if rep_type == 'upgma':
        file_path = out_path + '/JSON/' + AMR_gene + '_upgma.html'
        second_line = '\t\td3.json("' + AMR_gene + '_upgma.json"' + ')\n'

    elif rep_type == 'surrogates':
        file_path = out_path + '/JSON/' + AMR_gene + '_surrogates.html'
        second_line = '\t\td3.json("' + AMR_gene + '_surrogates.json"' + ')\n'

    else: # rep_type == 'standard'
        file_path = out_path + '/JSON/' + AMR_gene + '.html'
        second_line = '\t\td3.json("' + AMR_gene + '.json"' + ')\n'


    # Make HTML index file with appropriate JSON
    with open(file_path, 'w') as html_outfile, open(sample_data_path + '/index.html') as template:
        for line in template:
            html_outfile.write(line)

        html_outfile.write('\n')
        html_outfile.write(second_line)
        html_outfile.write('\t\t\t.then(data => {\n')
        html_outfile.write('\t\t\t\tdiv.selectAll("div")\n')
        html_outfile.write('\t\t\t\t\t.data([data])\n')
        html_outfile.write('\t\t\t\t\t.join("div")\n')
        html_outfile.write('\t\t\t\t\t.call(chart)\n\n')
        html_outfile.write('\t\t\t\tlet svg = div.select("svg")\n')
        html_outfile.write('\t\t\t\td3.select("#btn-save-svg")\n')
        html_outfile.write('\t\t\t\t\t.on("click", () => {\n')
        html_outfile.write('\t\t\t\t\t\tconst blob = serialise(svg)\n')
        html_outfile.write('\t\t\t\t\t\tdownload(blob, "clinker.svg")\n')
        html_outfile.write('\t\t\t\t\t})\n')
        html_outfile.write('\t\t\t})\n')
        html_outfile.write('\t</script>\n')
        html_outfile.write('</html>')


def get_cluster_data_genes_uid_list(json_cluster_data, genome_ids):
    """
    Given cluster data for a JSON clustermap representation of a gene's neighborhoods and a list of genomes to extract
    from, returns a list containing every gene uid present in the neighborhoods.
    """
    gene_uids_list = []
    for genome_id in genome_ids:
        for gene_key, gene_data in json_cluster_data[genome_id]["loci"]["genes"].items():
            gene_uids_list.append(json_cluster_data[genome_id]["loci"]["genes"][gene_key]["uid"])
    return gene_uids_list


def remove_defunct_clustermap_data(json_data):
    """
    For UPGMA representative clustermap JSON representations created based on original JSON file, after removing all
    cluster data for genomes no longer included, this function allows us to delete residual links and gene groups that
    no longer apply.
    """
    genome_ids = []
    for cluster in json_data["clusters"]:
        genome_ids.append(cluster["name"])
    print("Genomes in clustermap UPGMA rep: {g}".format(g=genome_ids))
    gene_links = []
    for link in json_data["links"]:

        # Get both genome names
        genome_contig_details = link["uid"].split('-')
        genome_1 = genome_contig_details[0].split('_')[0]
        genome_2 = genome_contig_details[1].split('_')[0]

        if genome_1 in genome_ids and genome_2 in genome_ids:
            gene_links.append(link)

    gene_uids_list = get_cluster_data_genes_uid_list(json_data["clusters"], genome_ids)
    gene_groups = []
    for group in json_data["groups"]:
        if group["uid"] in gene_uids_list:
            gene_groups.append(group)

    json_data["links"] = gene_links
    json_data["groups"] = gene_groups

    return json_data