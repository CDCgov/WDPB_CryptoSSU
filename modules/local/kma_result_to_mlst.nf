process KMA_RESULT_TO_MLST {

    tag { "$meta.id " }


    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::pandas conda-forge::biopython" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'quay.io/biocontainers/biopython' }"



    // executor 'local'

    //publishDir "${params.outdir}/${prefix}", pattern: "${prefix}_{cgmlst,locus_qc}.csv", mode: 'copy'

    input:
    tuple val(meta), path(kma_result)
    val(scheme)

    output:
    tuple val(meta), path("*_cgmlst.csv"), emit: mlst
    tuple val(meta), path("*_locus_qc.csv"), emit: mlst_qc
    
    script:
    def prefix=task.ext.prefix?:"${meta.id}"
    
    """
    
    #ln -s ${scheme}.name .

    KMADIR=\$(dirname \$(realpath ${scheme[0]}))
    
    kma_result_to_mlst.py \
      "${kma_result}" \
      --alleles \$KMADIR/kma_alleles.name \
      --sample-id "${prefix}" \
      --locus-allele-delimiter "_" \
      -o ${prefix}_cgmlst.csv \
      > ${prefix}_locus_qc.csv
    """
}

