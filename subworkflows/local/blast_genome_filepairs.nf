include { DIAMOND_MAKEDB } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTP } from '../../modules/nf-core/diamond/blastp/main'

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

        //Split CSV channels to get pairs of files for All-vs-All BLAST using DIAMOND_BLASTP
        blastPairs = filepairs_ch
            .splitCsv(header: true)
            .map{ row ->
                tuple(row.meta, row.fasta)
            }

        dbPairs = filepairs_ch
            .splitCsv(header: true)
            .map{ row ->
                row.db
            }

        // Run DIAMOND_BLASTP as BLAST All-vs-All for all combinations of fasta and db
        blastResults = DIAMOND_BLASTP(blastPairs, dbPairs, "txt", blast_columns).txt.collect()

    emit:
        db_files = dbFiles
        blast_files = blastResults

}
