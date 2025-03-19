process COUNT_CALLED_ALLELES {

    tag { "$meta.id "  }

    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::pandas conda-forge::biopython" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'quay.io/biocontainers/biopython' }"



    // executor 'local'

   // publishDir "${params.outdir}/${prefix}", pattern: "${prefix}_called_allele_count.csv", mode: 'copy'

    input:
    tuple val(meta), path(cgmlst)

    output:
    tuple val(meta), path("*_called_allele_count.csv"), emit: allele_count
    
    script:
    def prefix=task.ext.prefix?:"${meta.id}"

    """
    count_called_alleles.py ${cgmlst} > ${prefix}_called_allele_count.csv
    """
}
