process EXTRACTION {
    tag "extraction"
    
    input:
      path rgi_path
      path gbk_path
      path output_path
      value num_neighbors
      value percent_cutoff

    output:
      path "${output_path}/fasta", emit: fasta_path
      path "${output_path}/blast", emit: blast_path

    script:
    """
    extraction.py $rgi_path $gbk_path $output_path -n $num_neighbors -p $percent_cutoff
    """
}
