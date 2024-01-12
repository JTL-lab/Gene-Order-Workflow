process SUBSET_GENOMES {
    tag "subsetting"
    label 'process_low'

    publishDir "${params.outdir}"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jtl-lab-dev/bioinf-workflows/gene-order-workflow' :
        'jtllab/gene-order-workflow' }"

    input:
      path panaroo_output_path
      path output_path
      val percent_cutoff

    output:
      path "${params.outdir}/subsetted_genomes.txt", emit: subsetted_genomes

    // This script is bundled with the pipeline, in nf-core/geneorderanalysis/bin
    script:
    """
    subset_genomes.py $panaroo_output_path $output_path -p $percent_cutoff
    """
}
