#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

/*
========================================================================================
    DEFINE PROCESSES
========================================================================================
*/

process NEIGHBORHOOD_EXTRACTION {
    input:
      path rgi_path
      path gbk_path
      path output_path

    output:
      path "${output_path}/fasta", emit: fasta_path

    script:
    """
    extraction.py $rgi_path $gbk_path $output_path
    """
}


process NEIGHBORHOOD_CLUSTERING {
    input:
      path fasta_path
      path blast_path
      path output_path

    output:
      path "${output_path}/clustering", emit: cluster_path

    script:
    """
    clustering.py $fasta_path $blast_path $output_path
    """
}

workflow{

   // Initialize channels from provided paths
   def rgi_ch = Channel.fromPath(params.rgi_path)
   def gbk_ch = Channel.fromPath(params.gbk_path)
   def output_ch = Channel.fromPath(params.output_path)
   def blast_ch = Channel.fromPath(params.blast_path)

   // Extract neighborhoods of fixed window size
   NEIGHBORHOOD_EXTRACTION(rgi_ch, gbk_ch, output_ch)

   // Get fasta path from module output and initialize channel
   def fasta_ch = NEIGHBORHOOD_EXTRACTION.out.fasta_path

   // Cluster extracted neighborhoods
   NEIGHBORHOOD_CLUSTERING(fasta_ch, blast_ch, output_ch)

}
