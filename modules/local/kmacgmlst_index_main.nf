process KMACGMLST_INDEX {
    tag "$allele_seq"
    label 'process_medium'
        conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma%3A1.4.9--he4a0461_2' :
        'biocontainers/kma%3A1.4.9--he4a0461_2' }"

    input:
    path(allele_seq)
    
    output:
    path "*.{b,name}"    , emit: index // kma_Crypto_allele.comp.b, kma_Crypto_alleles.length.b, kma_Crypto_alleles.name, kma_Crypto_alleles.seq.b
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    
        kma index -i $allele_seq -o kma_alleles \\
        $args \\
        &> kmacgmlst_index.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmacgmlst_index: \$(kmacgmlst_indexversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
