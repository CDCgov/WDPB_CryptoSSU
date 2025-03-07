process SKESA {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::skesa=2.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/skesa:2.5.1--hdcf5f25_0' :
        'biocontainers/skesa:2.5.1--hdcf5f25_0' }"

    input:
    tuple val(meta), path(shortreads)

    output:
    tuple val(meta), path('*.fasta.gz')       , emit: contigs
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def short_reads = shortreads ? ( meta.single_end ? "$shortreads" : "${shortreads[0]},${shortreads[1]}" ) : ""
    """
    skesa \\
        --cores $task.cpus \\
        --reads $short_reads \\
        $args | \\
        gzip > ${prefix}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skesa: \$(echo \$(skesa --version 2>&1) | tail -n 1 | cut -f4 -d" ")
    END_VERSIONS
    """
}
