include { ALLELES_REFORMAT                      } from '../../modules/local/alleleDB_reformat'
include { REFERENCE_ALLELES                     } from '../../modules/local/reference_alleles'
include { ETOKI_DB_INDEX                        } from '../../modules/local/etoki_db_index'

workflow PREPARE_ALLELEDB{

    take:
    alleleDB_reformat     //channel: path to allele database

    main:
    ch_versions = Channel.empty()

    ALLELES_REFORMAT(alleleDB_reformat)
    ch_versions=ch_versions.mix(ALLELES_REFORMAT.out.versions)
    
    REFERENCE_ALLELES(alleleDB_reformat)
    ch_versions=ch_versions.mix(REFERENCE_ALLELES.out.versions)

    ETOKI_DB_INDEX(ALLELES_REFORMAT.out.alleles_reformat,REFERENCE_ALLELES.out.reference_alleles)
    ch_versions=ch_versions.mix(ETOKI_DB_INDEX.out.versions)


    emit:
    alleles_reformat    =     ALLELES_REFORMAT.out.alleles_reformat
    reference_alleles   =     REFERENCE_ALLELES.out.reference_alleles
    etoki_md5sum        =     ETOKI_DB_INDEX.out.etoki_md5sum

    versions            =     ch_versions


}
