//
// Check GBK filepairs samplesheet and get read channels in order to BLAST each pair using DIAMOND
//
include { DIAMOND_BLASTP } from '../../modules/nf-core/modules/nf-core/diamond/blastp/main'
include { DIAMOND_MAKEDB } from '../../modules/nf-core/modules/nf-core/diamond/makedb/main'

workflow BLAST_GENOME_FILEPAIRS {
    take:
        genomes
        genome_pairs

    main:
        // Specify BLAST output file columns
        def blast_columns = "qseqid sseqid pident bitscore"

        // Make BLAST DB for every unique genome assembly
        DIAMOND_MAKEDB(genomes)
        //DIAMOND_BLASTP(genome_paths_ch, blast_dbs, "txt", blast_columns).txt.set{ blast_file_paths }

    emit:
        blast_dbs = DIAMOND_MAKEDB.out.db
        //versions = BLAST_GENOME_FILEPAIRS.out.versions // channel: [ versions.yml ]
}
