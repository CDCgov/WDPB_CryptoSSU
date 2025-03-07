process QUAST {
    //tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321heaaa4ec_4' :
        'biocontainers/quast:5.2.0--py39pl5321heaaa4ec_4' }"

    input:
    //tuple val(meta) , path(consensus)
    //tuple val(meta2), path(fasta)
    //tuple val(meta3), path(gff)
    path(contigs)
    //path(contigs)
    path(fasta)
    path(gff)
    val use_fasta
    val use_gff

    output:
    path("${prefix}")           , emit: results
    path("*.tsv")               , emit: tsv
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    prefix        = task.ext.prefix ?: 'quast'
    def features  = gff             ?  "--features $gff" : ''
    def reference = fasta           ?  "-r $fasta"       : ''
    """
    quast.py \\
        --output-dir $prefix \\
        $reference \\
        $features \\
        --threads $task.cpus \\
        $args \\
        ${contigs.join(' ')}

    #ln -s ${prefix}/report.tsv ${prefix}.tsv
     ln -s ${prefix}/report.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """

}
