process CRYPTO_SSU_CLASSIFICATION_skesa {
    tag "$meta.id"
    label 'process_medium'




    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path('*.18S.txt'), emit: blast_results_skesa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #zcat $fasta > ${prefix}.fasta
    mkdir ./temp
    cp $fasta ./temp
    #gzip -d --force ./$fasta
    gunzip ./temp/$fasta
    mkdir ./${prefix}

    source /apps/x86_64/miniconda3/20230728/bin/activate /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/conda/py-blast

    18S_toolv06.py \\
        --query=./temp/${prefix}.fasta \\
        --reference_folder=${projectDir}/assets/databases/18S/ \\
        --resultsdir=./ \\
        --localdir=./${prefix}
        $args

    #mv ${prefix}/sorted_blastresults/blast_csv/${prefix}.csv ${prefix}.18S.txt
     mv ${prefix}.18SResults.csv ${prefix}.18S.txt
    """
}
