process CLUSTERING {
    tag "clustering"
    label 'process_high'
    label 'process_high_memory'

    publishDir "${params.output_path}"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity') {
        container "library://jtl-lab-dev/bioinf-workflows/gene-order-workflow"
    }

    input:
      path faa_path
      path fasta_path
      path blast_path
      path output_path
      val num_neighbors
      val inflation
      val epsilon
      val minpts

    output:
      path "${output_path}/clustering", emit: cluster_path

    script:
    """
    clustering.py $faa_path $fasta_path $blast_path $output_path -n $num_neighbors -i $inflation -e $epsilon -m $minpts
    """
}
