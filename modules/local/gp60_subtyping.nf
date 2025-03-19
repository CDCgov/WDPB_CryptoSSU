process SSU_CLASSIFICATION_GP60 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.32':
        'biocontainers/bioperl:v1.7.2-3-deb_cv1' }"
    input:
    tuple val(meta), path(blast_results)
    tuple val(meta), path(fasta)
    

    output:
    path('*.gp60.txt'), emit: gp60

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    gp60Typer_nf-core.pl \\
        --fasta ${prefix}.fasta \\
        --blastresultsgp60 $blast_results \\
        --data wgs \\
        --out ${prefix}.gp60.txt
        $args 
    """
}
