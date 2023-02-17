//
// Check GBK filepairs samplesheet and get read channels in order to BLAST each pair using DIAMOND
//
include { DIAMOND_BLASTP } from '../../modules/nf-core/modules/nf-core/diamond/blastp/main'
include { DIAMOND_MAKEDB } from '../../modules/nf-core/modules/nf-core/diamond/makedb/main'

workflow RUN_DIAMOND {
    take:
    csv_path

    main:
    Channel
        .fromPath( csv_path )
        .splitCsv( header: true, sep: ',')
        .map { row-> tuple(path(row.blast_dir), file(row.genome_1), file(row.genome_2)) }
        .DIAMOND_MAKEDB()
        .DIAMOND_BLASTP()

    emit:
    versions = RUN_DIAMOND.out.versions // channel: [ versions.yml ]
}

process MAKE_FILEPAIRS_DB {
    tag "$make_filepairs_db"
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
    tuple path(blast_dir), file(fasta_1), file(fasta_2)

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
