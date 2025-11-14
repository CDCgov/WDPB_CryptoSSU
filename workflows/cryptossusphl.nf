def valid_params =[assemblers : ['skesa','unicycler']]
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowCryptossusphl.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { FASTQ_TRIM_FASTP_FASTQC } from '../subworkflows/local/fastq_trim_fastp_fastqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { KRAKEN2                     } from '../modules/nf-core/kraken2/kraken2/main'
include { UNTAR as UNTAR_KRAKEN2_DB   } from '../modules/nf-core/untar/main'
include { KRAKEN2_BUILD               } from '../modules/local/kraken2_build'
include { UNICYCLER                   } from '../modules/nf-core/unicycler/main'
include { SKESA                       } from '../modules/local/skesa'
include { ETOKI_MLST                  } from '../modules/local/etoki_mlst'
include { BLAST_MAKEBLASTDB as MAKE_18S   } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_MAKEBLASTDB as MAKE_GP60  } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN as RUN_18S            } from '../modules/nf-core/blast/blastn/main'
include { BLAST_BLASTN as RUN_GP60           } from '../modules/nf-core/blast/blastn/main'
include { SSU_CLASSIFICATION_18S      } from '../modules/local/ssu_classification'
include { SSU_CLASSIFICATION_GP60     } from '../modules/local/gp60_subtyping'
include { CAT_18S                     } from '../modules/local/cat_18S'
include { CAT_GP60                    } from '../modules/local/cat_gp60'
include { QUAST                       } from '../modules/nf-core/quast/main'
include { EXTRACT_QUAST               } from '../modules/local/extract_quast'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CRYPTOSSUSPHL {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input),params.platform
    )
    .sample_info
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)



    //
    // MODULE: Run FastQC
    //

    FASTQC (
        ch_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // 
    // Module: Read QC and trim adapters
    FASTQ_TRIM_FASTP_FASTQC(
         ch_fastq,
         [],
         params.save_trimmed_fail,
         false
      )
      ch_trimmed_fastq = FASTQ_TRIM_FASTP_FASTQC.out.reads
      ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)

    // Module: Kraken reads

if (params.kraken2_db) {
            if (params.kraken2_db.endsWith('.tar.gz')) {
                UNTAR_KRAKEN2_DB (
                    [ [:], params.kraken2_db ]
                )
                ch_kraken2_db = UNTAR_KRAKEN2_DB.out.untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_KRAKEN2_DB.out.versions)
            } else {
                ch_kraken2_db = Channel.value(file(params.kraken2_db))
            }
        } else {
            KRAKEN2_BUILD (
                params.kraken2_db_name
            )
            ch_kraken2_db = KRAKEN2_BUILD.out.db.first()
            ch_versions   = ch_versions.mix(KRAKEN2_BUILD.out.versions)
        }


     ch_kraken2_multiqc = Channel.empty()

     KRAKEN2(ch_trimmed_fastq,[],[],[])
    ch_kraken2_multiqc = KRAKEN2.out.report
     ch_versions        = ch_versions.mix(KRAKEN2.out.versions.first().ifEmpty(null))
     //ch_kraken_quast = Channel.empty()


    // Module: Unicycler or Skesa
    ch_quast_contigs = Channel.empty()
    if ("${params.assemblers}" == 'unicycler'){
   
    UNICYCLER(ch_trimmed_fastq)
    ch_contigs = UNICYCLER.out.contigs
    ch_quast_contigs = UNICYCLER.out.contigs.collect{it[1]}
    } else { ("${params.assemblers}" == 'skesa')
    SKESA(ch_trimmed_fastq)
    ch_contigs = SKESA.out.contigs
    ch_versions = ch_versions.mix(SKESA.out.versions.first())
    ch_quast_contigs = SKESA.out.contigs.collect{it[1]}
    }

   // Module: Quast on assembly
  
   

   QUAST(ch_quast_contigs,params.fasta,params.gff,true,params.gff)
   ch_quast_report=Channel.empty()
   EXTRACT_QUAST(QUAST.out.tsv) //,"${params.assemblers}"
   ch_quast_report=EXTRACT_QUAST.out.extracted_quast
   ch_quast_report.view()


   // Modules: Blast 18S and GP60
   // make blast db
    MAKE_18S(params.db_18s)
    MAKE_GP60(params.db_gp60)
    
    ch_blast_results_18S = Channel.empty()
    ch_blast_results_18S = RUN_18S(ch_contigs,MAKE_18S.out.db)

    ch_ssu_results_18S = Channel.empty() 
    ch_ssu_results_18S = SSU_CLASSIFICATION_18S(RUN_18S.out.blast_results)
    ch_ssu_results_18S.view()
    ch_all_blast_results = ch_ssu_results_18S.map { string -> string }.collect()
    ch_all_blast_results.view()
    CAT_18S(ch_all_blast_results)

    ch_blast_results_gp60 = Channel.empty()
    ch_blast_results_gp60 = RUN_GP60(ch_contigs,MAKE_GP60.out.db)

    ch_ssu_results_gp60 = Channel.empty()
    ch_ssu_results_gp60 = SSU_CLASSIFICATION_GP60(RUN_GP60.out.blast_results,ch_contigs)
    ch_all_gp60_results = ch_ssu_results_gp60.map { string -> string }.collect()
    ch_all_gp60_results.view()
    CAT_GP60(ch_all_gp60_results)

   


    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCryptossusphl.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCryptossusphl.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.versions.first().ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CAT_GP60.out.gp60_results.collectFile(name: 'results_gp60_mqc.csv'))
    ch_multiqc_files = ch_multiqc_files.mix(CAT_18S.out.ssu18S_results.collectFile(name: 'results_18S_mqc.csv'))
    ch_multiqc_files = ch_multiqc_files.mix(EXTRACT_QUAST.out.extracted_quast.collectFile(name: 'extracted_quast_mqc.csv'))
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
