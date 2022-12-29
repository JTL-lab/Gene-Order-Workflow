//
// Check GBK filepairs samplesheet and get read channels in order to BLAST each pair using DIAMOND
//
include { BLASTP } from '../modules/nf-core/modules/diamond/blastp'
include { MAKE_GENOME_FILEPAIRS } from '../modules/local/make_genome_filepairs'

process MAKE_FILEPAIRS_DB {
    tag "$make_filepairs_db"
    debug true
    q
    // Dimaond is limited to v2.0.9 because there is not a
    // singularity version higher than this at the current time.
    conda (params.enable_conda ? 'bioconda::diamond=2.0.9' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0'
    } else {
        container 'quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0'
    }
        
    input:
    tuple path(blast_subdir), file(fasta_1), file(fasta_2)
    
    script:
    """
    diamond \\
        makedb \\
        --threads $task.cpus \\
        --in  $fasta_1 \\
        -d $fasta_1 \\
        $options.args
    """
}

process BLAST_GENOME_FILEPAIRS {
    tag "$blast_genome_filepairs"
    debug true
    
    // Dimaond is limited to v2.0.9 because there is not a
    // singularity version higher than this at the current time.
    conda (params.enable_conda ? 'bioconda::diamond=2.0.9' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0'
    } else {
        container 'quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0'
    }
    
    input:
    tuple val(meta), path(fasta)
    path db

    script:
    def outfmt = '6';
    def out_ext = 'txt';
    """
    DB=`find -L ./ -name "*.dmnd" | sed 's/.dmnd//'`

    diamond \\
        blastp \\
        --threads $task.cpus \\
        --db \$DB \\
        --query $fasta \\
        --outfmt ${outfmt} ${columns} \\
        $args \\
        --out ${prefix}.${out_ext}
    """
}

workflow (
    //def fasta_path = params.fasta_path
    //def blast_path = params.blast_path 
    //def output_path = params.output_path
    
    Channel.fromPath(params.fasta_path \
        | splitCsv(header:true) \
        | map { row-> tuple(path(row.blast_subdir), file(row.genome_1), file(row.genome_2)) } \
        | MAKE_FILEPAIRS_DB 
        | BLAST_GENOME_FILEPAIRS
)
