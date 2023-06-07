process CLUSTERING {
    tag "clustering"
    //label 'process_high'
    //label 'process_high_memory'

    publishDir "${params.outdir}"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jtl-lab-dev/bioinf-workflows/gene-order-workflow' :
        'jtllab/gene-order-workflow' }"

    input:
      tuple val(meta), path('*.txt')
      path faa_path
      path fasta_path
      path blast_path
      path output_path
      val num_neighbors
      val inflation
      val epsilon
      val minpts

    output:
      path "${params.outdir}/clustering", emit: cluster_path

    // This script is bundled with the pipeline, in nf-core/geneorderanalysis/bin
    script:
    """
    clustering.py $faa_path $fasta_path $blast_path $output_path -n $num_neighbors -i $inflation -e $epsilon -m $minpts
    """
}
