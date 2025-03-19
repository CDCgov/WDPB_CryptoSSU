process EXTRACT_QUAST {
    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::pandas conda-forge::biopython" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas' }"
    
    input:
    path(quastreport)
    //val assemblers 

    output:
    path "*.csv"        , emit: extracted_quast
    path "versions.yml"     , emit:versions

    script:
    // def args = task.ext.args   ?: ''
    def report = "${projectDir}/quast_assembly/quast/report.tsv"

    """
    quast_extract.py \\
    $quastreport \\
    ./extracted_quast_mqc.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
