process SSU_CLASSIFICATION_18S {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3 conda-forge::pandas conda-forge::biopython" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas' }"


    input:
    tuple val(meta), path(blast_results)

    output:
    path('*.txt'), emit: results_ssu_18S_report

    when:
    task.ext.when == null || task.ext.when

   script:
    def args = task.ext.args ?: ''
    def inputdir = "${meta.id}"
    def tempdir = "temp_csv"
    def prefix = task.ext.prefix ?: "${meta.id}"
    
   """
    mkdir $inputdir
    cp $blast_results $inputdir
    python3 ${projectDir}/bin/18S_toolv1.2.py \\
        --inputdir $inputdir \\
        --resultsdir ./ \\
        --localdir $tempdir
     
    #mv ${prefix}/sorted_blastresults/blast_csv/${prefix}.csv ${prefix}.18S.txt
    mv ${prefix}.18SResults.csv ${prefix}.18S.txt
    
  """  
}
