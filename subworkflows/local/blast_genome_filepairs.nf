// Import the modules
include { DIAMOND_MAKEDB } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTP } from '../../modules/nf-core/diamond/blastp/main'

// Define the subworkflow
workflow BLAST_GENOME_FILEPAIRS {

    take:
        csv_ch
        assemblyFiles

    main:
        // Keep only columns required for clustering module
        def blast_columns = "qseqid sseqid pident bitscore"

        // Split CSV channel to get filepaths to assemblies
        csv_ch
            .splitCsv(header: true)
            .map { it.file_path }
            .set{ genomeFiles }

        // Create a database for each genome using nf-core MAKEDB module
        DIAMOND_MAKEDB(genomeFiles).db.set { dbFiles }

        // BLAST every assembly file against every DB for All-vs-All BLAST
        assemblyFiles.each { fasta_tuple ->
            DIAMOND_BLASTP(fasta_tuple, dbFiles, "txt", blast_columns)
        }

    emit:
        blastResults = DIAMOND_BLASTP.out.txt
}