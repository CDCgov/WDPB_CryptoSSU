process CAT_18S_skesa {
    label 'process_low'

    conda "conda-forge::perl=5.32.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    path(files_in)

    output:
    path("*_18S_mqc.csv"), emit: ssu18S_results_skesa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_list = files_in.collect { it.toString() }
    """
    echo -e "Sample Name,query_genome,db_bestmatch,species,pident,alignment_length,coverage,query_start,query_end,subject_start,subject_end,bitscore,NCE" > results_18S_mqc.csv
    cat ${file_list.join(' ')} | grep -v "^sample_name" | sed "s/\\t/,/g" >> results_18S_mqc.csv
    """
}
