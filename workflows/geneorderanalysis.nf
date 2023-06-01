/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowGeneorderanalysis.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input_file_path, params.assembly_path, params.gbk_path, params.output_path ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Consisting of a mix of local and nf-core modules, plus subworkflow for BLAST using nf-core/diamond
//
include { EXTRACTION } from '../modules/local/extraction'
include { MAKE_GENOME_SAMPLESHEET } from '../modules/local/make_genome_samplesheet'
include { MAKE_GENOME_FILEPAIRS } from '../modules/local/make_genome_filepairs'
include { CLUSTERING } from '../modules/local/clustering'
include { DIAMOND_BLASTP } from '../modules/nf-core/diamond/blastp/main'
include { DIAMOND_MAKEDB } from '../modules/nf-core/diamond/makedb/main' addParams( options: [args:'--task blastp --outfmt 6 --max-hsps 1'] )
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { BLAST_GENOME_FILEPAIRS } from '../subworkflows/local/blast_genome_filepairs'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENEORDERANALYSIS {

    // Initialize channels from provided paths
    def input_file_ch = Channel.fromPath(params.input_file_path)
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
        input_file_ch,
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

    def csv_ch = MAKE_GENOME_SAMPLESHEET.out.genome_paths_csv

    MAKE_GENOME_FILEPAIRS (
        assembly_ch,
        output_ch
    )

    def filepairs_ch = MAKE_GENOME_FILEPAIRS.out.genome_filepairs_csv

    //
    // MODULE: Run nf-core/diamond to obtain BLAST results
    //
    //def assemblyFiles = Channel.fromPath("${params.assembly_path}/*.{fa,faa,fna}")

    BLAST_GENOME_FILEPAIRS(csv_ch, filepairs_ch)

    def blastResults_ch = BLAST_GENOME_FILEPAIRS.out.blastResults

    //
    // MODULE: Run clustering
    //
    CLUSTERING (
        blastResults_ch,
        assembly_ch,
    	EXTRACTION.out.fasta_path,
    	EXTRACTION.out.blast_path,
    	output_ch,
    	num_neighbors,
    	inflation,
    	epsilon,
    	minpts
    )

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
