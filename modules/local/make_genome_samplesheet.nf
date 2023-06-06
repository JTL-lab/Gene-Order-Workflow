process MAKE_GENOME_SAMPLESHEET {
    tag "$make_genome_samplesheet"
    label 'process_low'

    publishDir "${params.output_path}"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jtl-lab-dev/bioinf-workflows/gene-order-workflow' :
        'jtllab/gene-order-workflow' }"

    input:
    path assembly_path
    path output_path

    output:
    path "${output_path}/genome_paths.csv", emit: genome_paths_csv

    // This script is bundled with the pipeline, in nf-core/geneorderanalysis/bin
    script:
    """
    make_genome_samplesheet.py $assembly_path $output_path
    """
}
