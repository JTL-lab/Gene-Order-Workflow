process MAKE_GENOME_FILEPAIRS {
    tag "$make_genome_filepairs"
    label 'process:low'

    publishDir "${params.output_path}"

    input:
    path assembly_path
    path output_path

    output:
    path "${output_path}/genome_pairs.csv", emit: csv_path

    script: // This script is bundled with the pipeline, in nf-core/geneorderanalysis/bin/
    """
    make_genome_filepairs.py $assembly_path $output_path
    """
}
