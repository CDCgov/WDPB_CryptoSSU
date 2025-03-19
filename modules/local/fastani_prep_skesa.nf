process FASTANI_PREP_SKESA {
    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::pandas conda-forge::biopython" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'quay.io/biocontainers/biopython' }"

    input:
    path(scaffolds)
    output:
    path '*skesa.txt', emit: fastani_prepfile_skesa
    path "versions.yml", emit:versions

    when: 
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def SCAFFOLDDIR_2 = "./assembly"
    
    """
    uDir=`dirname $SCAFFOLDDIR_2`
    

   prep_fastani_v2_skesa.py \\
	-assemblies \$uDir \\
	-output .
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
	


