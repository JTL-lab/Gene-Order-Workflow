process CLUSTERING {
    tag "clustering"
    input:
      path faa_path
      path fasta_path
      path blast_path
      path output_path
      val num_neighbors

    output:
      path "${output_path}/clustering", emit: cluster_path

    script:
    """
    clustering.py $faa_path $fasta_path $blast_path $output_path -n $num_neighbors
    """
}
