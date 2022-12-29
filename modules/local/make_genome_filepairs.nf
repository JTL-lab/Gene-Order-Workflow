process MAKE_GENOME_FILEPAIRS {
    tag "$make_genome_filepairs"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path fasta_path
    path blast_path
    path output_path

    output:
    path "${output_path}/genome_pairs.csv", emit: csv_path

    script: // This script is bundled with the pipeline, in nf-core/geneorderanalysis/bin/
    """
    make_genome_filepairs.py $fasta_path $blast_path $output_path
    """
}
