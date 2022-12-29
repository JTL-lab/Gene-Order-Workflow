/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowGeneorderanalysis.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.rgi_path, params.gbk_path, params.output_path ]
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
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { EXTRACTION } from '../modules/local/extraction'
include { MAKE_GENOME_FILEPAIRS } from '../modules/local/make_genome_filepairs'
include { RUN_DIAMOND } from '../subworkflows/local/run_diamond'
include { CLUSTERING } from '../modules/local/clustering'

include { DIAMOND_BLASTP } from '../modules/nf-core/modules/diamond/blastp/main'  addParams( options: [args:'--task blastp --outfmt 6 --max-hsps 1'] )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow GENEORDERANALYSIS {

    //ch_versions = Channel.empty()

    // Initialize channels from provided paths
    def rgi_ch = Channel.fromPath(params.rgi_path)
    def gbk_ch = Channel.fromPath(params.gbk_path)
    def output_ch = Channel.fromPath(params.output_path)

    //
    // MODULE: Run extraction 
    //
    EXTRACTION (
    	rgi_ch,
    	gbk_ch,
    	output_ch
    )
    
    //
    // MODULE: Run nf-core/diamond to obtain BLAST results
    //
    def fasta_path = EXTRACTION.out.fasta_path
    def blast_path = EXTRACTION.out.blast_path
    
    MAKE_GENOME_FILEPAIRS (
        fasta_path,
        blast_path,
        output_ch
    )
    
    csv_path = MAKE_GENOME_FILEPAIRS.out.csv_path
    
    //def fasta_filepaths = Channel.from('$fasta_path/*.fasta', chec)
    //RUN_DIAMOND (
    //    csv_path,
    //	fasta_path,
    //	blast_path
    //)

    //DIAMOND_BLASTP (
    //	tuple(val/path),
    //	path db,
    //	val out_ext
    //	val blast_columns
    //)

    //
    // MODULE: Run clustering
    //
    
    CLUSTERING (
    	EXTRACTION.out.fasta_path,
    	EXTRACTION.out.blast_path,
    	output_ch
    )
    
    VISUALIZATION (
    
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

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
