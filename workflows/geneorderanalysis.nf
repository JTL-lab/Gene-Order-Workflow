/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowGeneorderanalysis.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.assembly_path, params.extract_path, params.gbk_path, params.output_path ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
//if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { EXTRACTION } from '../modules/local/extraction'
include { MAKE_GENOME_SAMPLESHEET } from '../modules/local/make_genome_samplesheet'
//include { MAKE_GENOME_FILEPAIRS } from '../modules/local/make_genome_filepairs'
//include { CREATE_DB } from '../modules/local/all_vs_all_blast.nf'
//include { ALL_VS_ALL_BLAST } from '../modules/local/all_vs_all_blast.nf'
include { CLUSTERING } from '../modules/local/clustering'
include { DIAMOND_BLASTP } from '../modules/nf-core/diamond/blastp/main'
include { DIAMOND_MAKEDB } from '../modules/nf-core/diamond/makedb/main' addParams( options: [args:'--task blastp --outfmt 6 --max-hsps 1'] )
include { INPUT_CHECK } from '../subworkflows/local/input_check'
//include { INPUT_GENOMES_CHECK } from '../subworkflows/local/input_genomes_check.nf'
//include { INPUT_GENOME_PAIRS_CHECK } from '../subworkflows/local/input_genome_pairs_check.nf'
include { BLAST_GENOME_FILEPAIRS } from '../subworkflows/local/blast_genome_filepairs'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
//def multiqc_report = []

workflow GENEORDERANALYSIS {

    //ch_versions = Channel.empty()

    // Initialize channels from provided paths
    def assembly_ch = Channel.fromPath(params.assembly_path)
    def extract_ch = Channel.fromPath(params.extract_path)
    def gbk_ch = Channel.fromPath(params.gbk_path)
    def output_ch = Channel.fromPath(params.output_path)

    // Optional extraction params
    num_neighbors = params.num_neighbors
    percent_cutoff = params.percent_cutoff

    // Optional clustering params
    inflation = params.inflation
    epsilon = params.epsilon
    minpts = params.minpts

    //
    // MODULE: Run extraction 
    //
    EXTRACTION (
        extract_ch,
    	gbk_ch,
    	output_ch,
    	num_neighbors,
    	percent_cutoff
    )

    //
    // MODULE: Make CSV listing path to every unique genome assembly for DIAMOND MAKEDB
    //
    MAKE_GENOME_SAMPLESHEET (
        assembly_ch,
        output_ch
    )

    csv_ch = MAKE_GENOME_SAMPLESHEET.out.genome_paths_csv

    //
    // MODULE: Make CSV listing all unique genome pairs for All-vs-All DIAMOND BLASTP
    //
    //MAKE_GENOME_FILEPAIRS (
    //    params.assembly_path,
    //    params.output_path
    //)

    //INPUT_GENOME_PAIRS_CHECK(MAKE_GENOME_FILEPAIRS.out.genome_filepairs_csv)

    //
    // MODULE: Run nf-core/diamond to obtain BLAST results
    //
    //BLAST_GENOME_FILEPAIRS(assembly_ch)
    //BLAST_GENOME_FILEPAIRS(INPUT_GENOMES_CHECK.out.genomes, INPUT_GENOME_PAIRS_CHECK.out.genome_pairs)
    def assemblyFiles = Channel.fromPath("${params.assembly_path}/*.{fa,faa,fna}")
    assemblyFiles.view()
    def modifiedChannel = assemblyFiles.collect { path ->
        tuple(path, path)
    }

    // Print the channel of tuples
    modifiedChannel.view()

    BLAST_GENOME_FILEPAIRS(csv_ch, modifiedChannel)

    // Create the databases for each genome
    //CREATE_DB(assembly_ch)

    //dbFiles = CREATE_DB.out.dbFiles

    // Perform All-vs-All BLAST for each genome
    //ALL_VS_ALL_BLAST(assembly_ch, dbFiles)

    //
    // MODULE: Run clustering
    //
    CLUSTERING (
        BLAST_GENOME_FILEPAIRS.out.blastResults,
        assembly_ch,
    	EXTRACTION.out.fasta_path,
    	EXTRACTION.out.blast_path,
    	output_ch,
    	num_neighbors,
    	inflation,
    	epsilon,
    	minpts
    )
    
    //workflow_summary    = WorkflowGeneorderanalysis.paramsSummaryMultiqc(workflow, summary_params)
    //ch_workflow_summary = Channel.value(workflow_summary)

    //ch_multiqc_files = Channel.empty()
    //ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    //MULTIQC (
    //    ch_multiqc_files.collect()
    //)
    //multiqc_report = MULTIQC.out.report.toList()
    //ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//workflow.onComplete {
//    if (params.email || params.email_on_fail) {
//        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//    }
//    NfcoreTemplate.summary(workflow, params, log)
//}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
