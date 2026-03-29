/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.input) {
    ch_input = file(params.input)
} else {
    exit 1, 'Input samplesheet not specified!'
}

if (!params.fasta) { exit 1, "CircRNA analysis requires '--fasta' (reference genome)." }
if (!params.gtf)   { exit 1, "CircRNA analysis requires '--gtf' (gene annotation)." }

if ((params.run_isocirc || params.run_cirilong) && !params.circrna_db) {
    exit 1, "isocirc and ciri-long require '--circrna_db' (circRNA database)."
}

if (!params.run_isocirc && !params.run_circfl && !params.run_cirilong && !params.run_circnick) {
    exit 1, "At least one tool must be active. Use --run_isocirc, --run_circfl, --run_cirilong, or --run_circnick."
}

if (params.run_circnick) {
    if (!params.circnick_species) {
        exit 1, "Parameter '--circnick_species' is required when '--run_circnick' is set. Valid options: 'mouse', 'human'"
    }
    if (params.circnick_species != 'mouse' && params.circnick_species != 'human') {
        exit 1, "Invalid --circnick_species: '${params.circnick_species}'. Valid options: 'mouse', 'human'"
    }
    if (!params.circnick_liftover_chain) {
        log.warn "No --circnick_liftover_chain provided. circnick-lrs will use built-in ${params.circnick_species == 'mouse' ? 'mm10' : 'hg19'} coordinates."
    }
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = params.multiqc_config \
    ? Channel.fromPath(params.multiqc_config, checkIfExists: true) \
    : Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: false)
ch_multiqc_custom_config = Channel.empty()

////////////////////////////////////////////////////
/* --    IMPORT MODULES/SUBWORKFLOWS           -- */
////////////////////////////////////////////////////

include { NANOLYSE            } from '../modules/nf-core/nanolyse/main'
include { FASTQC              } from '../modules/nf-core/fastqc/main'
include { NANOPLOT            } from '../modules/nf-core/nanoplot/main'
include { MULTIQC             } from '../modules/nf-core/multiqc/main'
include { FASTQ_RENAME        } from '../modules/local/fastq_rename'
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { CIRCRNA_ANALYSIS    } from '../subworkflows/local/circrna_analysis'

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow NANOCIRC {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( ch_input )

    // Ensure all input files are named *.fastq.gz — required by nf-core NanoPlot
    // module which detects input type from extension. Input is already validated
    // as .gz by INPUT_CHECK; this just standardises the filename.
    FASTQ_RENAME ( INPUT_CHECK.out.fastq )
    ch_fastq = FASTQ_RENAME.out.fastq

    /*
     * MODULE: DNA contaminant removal using NanoLyse (optional)
     */
    if (params.run_nanolyse) {
        if (!params.nanolyse_db) {
            exit 1, "Parameter '--nanolyse_db' is required when '--run_nanolyse' is set."
        }
        NANOLYSE ( ch_fastq, file(params.nanolyse_db, checkIfExists: true) )
        ch_fastq             = NANOLYSE.out.fastq
        ch_software_versions = ch_software_versions.mix(NANOLYSE.out.versions.first().ifEmpty(null))
    }

    /*
     * QC: FastQC and NanoPlot — both handle .fastq.gz natively
     */
    ch_fastqc_multiqc = Channel.empty()
    if (!params.skip_qc) {
        if (!params.skip_fastqc) {
            FASTQC ( ch_fastq )
            ch_software_versions = ch_software_versions.mix(FASTQC.out.versions.first().ifEmpty(null))
            ch_fastqc_multiqc    = FASTQC.out.zip.collect { it[1] }.ifEmpty([])
        }
        if (!params.skip_nanoplot) {
            NANOPLOT ( ch_fastq )
            ch_software_versions = ch_software_versions.mix(NANOPLOT.out.versions.first().ifEmpty(null))
        }
    }

    /*
     * SUBWORKFLOW: circRNA detection, BED12 conversion, merge and confidence scoring
     * Tools decompress inline if they require unzipped input.
     */
    CIRCRNA_ANALYSIS (
        ch_fastq,
        file(params.fasta,       checkIfExists: true),
        file(params.gtf,         checkIfExists: true),
        params.circrna_db ? file(params.circrna_db, checkIfExists: true) : file('NO_FILE')
    )
    ch_software_versions = ch_software_versions.mix(CIRCRNA_ANALYSIS.out.versions.ifEmpty(null))

    /*
     * MODULE: MultiQC
     */
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowNanocirc.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_fastqc_multiqc.mix(
                ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
            ).collect().ifEmpty([]),
            ch_multiqc_config.collect().ifEmpty([]),
            ch_multiqc_custom_config.collect().ifEmpty([]),
            []
        )
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions)
    }
}

////////////////////////////////////////////////////
/* --           COMPLETION EMAIL               -- */
////////////////////////////////////////////////////

workflow.onComplete {
    NfcoreTemplate.summary(workflow, params, log)
    if (params.email) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, [])
    }
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
