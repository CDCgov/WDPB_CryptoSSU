process ETOKI_DB_INDEX {
    // tag "$fasta"
    label 'process_medium'



    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/etoki%3A1.2.3--hdfd78af_0' :
        'biocontainers/fastani:etoki%3A1.2.3--hdfd78af_0' }"


    input:
    path(alleles_reformat)
    path(reference_alleles)
    
    output:
    path "*.csv"            , emit: etoki_md5sum 
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
       
        EToKi.py MLSTdb \\
        -i $alleles_reformat \\
        -r $reference_alleles \\
        -d etoki_index_md5hash.csv \\
        $args \\
        &> etoki_index.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        etoki_index: \$(etoki_index_indexversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
