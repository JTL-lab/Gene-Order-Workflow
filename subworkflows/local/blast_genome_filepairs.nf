// Import the modules
include { DIAMOND_MAKEDB } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTP } from '../../modules/nf-core/diamond/blastp/main'

// Define the subworkflow
workflow BLAST_GENOME_FILEPAIRS {

    take:
        csv_ch
        filepairs_ch

    main:
        // Keep only columns required for clustering module
        def blast_columns = "qseqid sseqid pident bitscore"

        // Split CSV channel to get filepaths to assemblies
        genomeFiles = csv_ch
            .splitCsv(header: true)
            .map { it.file_path }

        // Create a database for each genome using nf-core MAKEDB module
        dbFiles = DIAMOND_MAKEDB(genomeFiles).db.collect()

        //Split filepairs_ch to get pairs of files for All-vs-All BLAST using DIAMOND_BLASTP
        blastPairs = filepairs_ch
            .splitCsv(header: true)
            .map{ row ->
                tuple(row.meta, row.fasta)
            }
            .view()

        dbPairs = filepairs_ch
            .splitCsv(header: true)
            .map{ row ->
                row.db
            }
            .view()

        // BLAST every assembly file against every DB for All-vs-All BLAST
        blastResults_ch = Channel.empty()

        blastPairs.each { fasta_tuple ->
            dbPairs.each { db_file ->
                DIAMOND_BLASTP(fasta_tuple, db_file, "txt", blast_columns).txt
                    .collect { blastResult ->
                        blastResults_ch << blastResult
                    }
            }
        }

        blastResults_ch.set { blastResults }

    emit:
        blastResults
}
