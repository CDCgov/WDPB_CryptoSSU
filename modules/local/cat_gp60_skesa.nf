process CAT_GP60_skesa {
    label 'process_low'

    conda "conda-forge::perl=5.32.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    path(files_in)

    output:
    path("*_gp60_mqc.tsv"), emit: gp60_results_skesa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_list = files_in.collect { it.toString() }
    """
    echo -e "Sample Name\tgp60_Subtype\tFasta_Header\tDatabase_Best_Match\tLength(bp)\tBlast_Identity\tCoverage\tNCE" > results_gp60_mqc.tsv
    awk -F'\t' 'BEGIN {OFS="\t"} FNR == 1 && NR != 1 { next } \$1 != "Sample Name" { \$8 = (\$8 == "" ? "NA" : \$8); print }' ${file_list.join(' ')} >> results_gp60_mqc.tsv
    """
}
