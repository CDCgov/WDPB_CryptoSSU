process CRYPTO_GP60_SUBTYPING_skesa {
    tag "$meta.id"
    label 'process_low'

    container "${baseDir}/assets/gp60.sif"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path('*.gp60.txt'), emit: gp60_skesa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dbname = "${db.simpleName}"
    """
    mkdir database
    cp \$(realpath $db).n* database/
    zcat $fasta > ${prefix}.fasta

    gp60Typer.pl \\
        --fasta ${prefix}.fasta \\
        --blastdb database/$dbname \\
        --out ${prefix}.gp60.txt \\
        $args 
    """
}
