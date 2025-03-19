process KMA_ALIGN {

    tag { "$meta.id "}
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma%3A1.4.9--he4a0461_2' :
        'biocontainers/kma%3A1.4.9--he4a0461_2' }"


    //publishDir "${params.outdir}/${prefix}", pattern: "${prefix}_kma*.{c,t}sv", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    val(scheme)

    output:
    tuple val(meta), path("*_kma.csv"), emit: res
    tuple val(meta), path("*_kma_mapstat.tsv"), emit: mapstat
    tuple val(meta), path("*_kma_align_provenance.yml"), emit: provenance

    script:
    def prefix=task.ext.prefix?:"${meta.id}"
    
    """

    printf -- "- process_name: kma_align\\n"       >> ${prefix}_kma_align_provenance.yml
    printf -- "  tools:\\n"                        >> ${prefix}_kma_align_provenance.yml
    printf -- "    - tool_name: kma\\n"            >> ${prefix}_kma_align_provenance.yml
    printf -- "      tool_version: \$(kma -v 2>&1 | cut -d '-' -f 2)\\n" >> ${prefix}_kma_align_provenance.yml
    printf -- "      parameters:\\n"               >> ${prefix}_kma_align_provenance.yml
    printf -- "        - parameter: -ef\\n"        >> ${prefix}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${prefix}_kma_align_provenance.yml
    printf -- "        - parameter: -cge\\n"       >> ${prefix}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${prefix}_kma_align_provenance.yml
    printf -- "        - parameter: -boot\\n"      >> ${prefix}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${prefix}_kma_align_provenance.yml
    printf -- "        - parameter: -1t1\\n"       >> ${prefix}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${prefix}_kma_align_provenance.yml
    printf -- "        - parameter: -mem_mode\\n"  >> ${prefix}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${prefix}_kma_align_provenance.yml
    printf -- "        - parameter: -and\\n"       >> ${prefix}_kma_align_provenance.yml
    printf -- "          value: null\\n"           >> ${prefix}_kma_align_provenance.yml
    
    #ln -s -f ${scheme}.comp.b .
    # ln -s -f ${scheme}.length.b .
    # ln -s -f ${scheme}.name .
    # ln -s -f ${scheme}.seq.b .
    mkdir ./${prefix}_kma 
    KMADIR=\$(dirname \$(realpath ${scheme[0]}))
    
    kma -o ./${prefix}_kma -t_db \$KMADIR/kma_alleles -ipe ${reads[0]} ${reads[1]} -ef -cge -mem_mode -t ${task.cpus} -1t1 -and -boot

    head -n 1 ${prefix}_kma.res | tr -d '#' | awk '{print tolower(\$0)}' | tr \$'\\t' ',' > ${prefix}_kma.csv
    tail -qn+2 ${prefix}_kma.res | tr -d ' ' | tr \$'\\t' ',' >> ${prefix}_kma.csv

    mv ${prefix}_kma.mapstat ${prefix}_kma_mapstat.tsv
    """
}
