process CAT_GP60 {
    label 'process_low'

    conda "conda-forge::perl=5.32.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    path(files_in)

    output:
    path("*_gp60_mqc.csv"), emit: gp60_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_list = files_in.collect { it.toString() }
    """
    #echo $files_in
    echo -e "Sample Name,gp60_Subtype,Fasta_Header,Database_Best_Match,Length(bp),Blast_Identity,Coverage,NCE" > results_gp60_mqc.csv
    cat ${file_list.join(' ')} | grep -v "^sample_name" | sed "s/\\t/,/g" >> results_gp60_mqc.csv
    #awk -F'\t' 'BEGIN {OFS="\t"} FNR == 1 && NR != 1 { next } \$1 != "Sample Name" { \$8 = (\$8 == "" ? "NA" : \$8); print }' ${file_list.join(' ')} >> results_gp60_mqc.csv
    """
}
