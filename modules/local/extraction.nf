process EXTRACTION {
    tag "extraction"
    //label 'process_medium'

    publishDir "${params.outdir}"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' ?
        'library://jtl-lab-dev/bioinf-workflows/gene-order-workflow' :
        'jtllab/gene-order-workflow' }"

    input:
      path input_file_path
      path extract_path
      path gbk_path
      path output_path
      val num_neighbors
      val percent_cutoff

    output:
      path "${params.outdir}/fasta", emit: fasta_path
      path "${params.outdir}/diamond", emit: blast_path

    // This script is bundled with the pipeline, in nf-core/geneorderanalysis/bin
    script:
    """
    extraction.py $input_file_path $extract_path $gbk_path $output_path -n $num_neighbors -p $percent_cutoff
    """
}
