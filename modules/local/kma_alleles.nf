process KMA_ALLELES {
    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::pandas conda-forge::biopython" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'quay.io/biocontainers/biopython' }"

    input:
    path(fasta)

    output:
    path '*.fasta'       , emit: fasta
    path "versions.yml", emit:versions

    //def args = task.ext.args   ?: ''
    """
    kma_alleles_format.py \\
	-bn $fasta \\
	-output .
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
	
