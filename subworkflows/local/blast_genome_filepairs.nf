//
// Check GBK filepairs samplesheet and get read channels in order to BLAST each pair using DIAMOND
//
include { DIAMOND_BLASTP } from '../../modules/nf-core/modules/nf-core/diamond/blastp/main'
include { DIAMOND_MAKEDB } from '../../modules/nf-core/modules/nf-core/diamond/makedb/main'

workflow BLAST_GENOME_FILEPAIRS{
    take:
        genome_paths_csv
        genome_filepairs_csv

    main:
        // Initialize channels to hold outputs
        blast_dbs = Channel.empty()
        blast_outputs = Channel.empty()

        // Specify BLAST output file columns
        def blast_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

        // Create channel to hold paths to all whole genome assembly files as values
        Channel
            .fromPath ( genome_paths_csv )
            .splitCsv ( header: true, sep: ',' )
            .map { row -> path(row.file_path) }
            .set { genome_paths_ch }

        // Create channel to hold paths to files to blast against each other
        Channel
            .fromPath( genome_filepairs_csv )
            .splitCsv( header: true, sep: ',' )
            .map { row-> tuple(path(row.blast_dir), file(row.genome_1), file(row.genome_2)) }
            .set { genome_pairs_ch }

        // Make BLAST DB for every unique genome assembly
        DIAMOND_MAKEDB(genome_paths_ch).db.set{ blast_dbs }
        DIAMOND_BLASTP(genome_paths_ch, blast_dbs, "txt", blast_columns).txt.set{ blast_file_paths }

    emit:
        blast_file_paths
        //versions = BLAST_GENOME_FILEPAIRS.out.versions // channel: [ versions.yml ]
}
