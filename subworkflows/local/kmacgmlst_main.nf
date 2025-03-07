#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { hash_files 	}                from '../../modules/local/hash_files.nf'
include { fastp         }                from '../../modules/local/fastp.nf'
include { KMA_ALIGN     }                from '../../modules/local/kma_align.nf'
include { KMA_RESULT_TO_MLST    }        from '../../modules/local/kma_result_to_mlst.nf'
include { COUNT_CALLED_ALLELES }         from '../../modules/local/count_called_alleles.nf'
// include { pipeline_provenance }          from './modules/local/provenance.nf'
// include { collect_provenance }           from './modules/local/provenance.nf'

if (params.profile){
	println("Profile should have a single dash: -profile")
	System.exit(1)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow KMACGMLST_CALL{
	take:
    ch_bbmap_trimmedreads
	ch_scheme
	// ch_workflow_metadata = Channel.value([
	// 	workflow.sessionId,
	// 	workflow.runName,
	// 	workflow.manifest.name,
	// 	workflow.manifest.version,
	// 	workflow.start,
	// ])
	// if (params.input != 'NO_FILE') {
	// 	ch_fastq = Channel.fromPath(params.ch_samplesheet).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
	// } else {
	// 	ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
	// }


	main:

	// hash_files(ch_fastq.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq-input")))
	// fastp(ch_fastq)
	KMA_ALIGN(ch_bbmap_trimmedreads, ch_scheme)
	KMA_RESULT_TO_MLST(KMA_ALIGN.out.res, ch_scheme)
	COUNT_CALLED_ALLELES(KMA_RESULT_TO_MLST.out.mlst)
        ch_versions = Channel.empty()

	// if (params.collect_outputs) {
	// fastp.out.csv.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_fastp.csv", storeDir: params.outdir, keepHeader: true, sort: { it -> it.readLines()[1].split(',')[0] })

	// count_called_alleles.out.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_called_allele_count.csv", storeDir: params.outdir, keepHeader: true, sort: { it -> it.readLines()[1].split(',')[0] })

	// kma_result_to_mlst.out.mlst.map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_cgmlst.csv", storeDir: params.outdir, keepHeader: true, sort: { it -> it.readLines()[1].split(',')[0] })
	// }

	// Collect Provenance
	// The basic idea is to build up a channel with the following structure:
	// [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
	// At each step, we add another provenance file to the list using the << operator...
	// ...and then concatenate them all together in the 'collect_provenance' process.
	// ch_sample_ids = ch_fastq.map{ it -> it[0] }
	// ch_provenance = ch_sample_ids
	// ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)
	// ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
	// ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
	// ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
	// ch_provenance = ch_provenance.join(kma_align.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

	// collect_provenance(ch_provenance)
	// emit:
    KMA_ALIGN.out.res 								// [val(meta), res]
    KMA_RESULT_TO_MLST.out.mlst  					// channel: [val(meta), "*_cgmlst.csv"]
    KMA_RESULT_TO_MLST.out.mlst_qc 					// channel: [ val(meta), "*_locus_qc.csv")]
	COUNT_CALLED_ALLELES.out.allele_count 			// channel: [ val(meta), "*_called_allele_count.csv"]
	versions = ch_versions   
          
}
