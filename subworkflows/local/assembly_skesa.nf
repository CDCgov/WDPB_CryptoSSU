// Assembly and downstream processing for Unicycler scaffolds
// Credit to: nf-core/viralrecon (https://github.com/nf-core/viralrecon)

include { SKESA                      } from '../../modules/local/skesa'
include { GUNZIP as GUNZIP_SCAFFOLDS } from '../../modules/nf-core/gunzip/main'
include { ASSEMBLY_QC_SKESA                } from './assembly_qc_skesa'

workflow ASSEMBLY_SKESA {
    take:
    reads        // channel: [ val(meta), [ reads ] ]
    fasta        // channel: /path/to/genome.fasta
    gff          // channel: /path/to/genome.gff
    blast_db     // channel: /path/to/blast_db/
    blast_header // channel: /path/to/blast_header.txt

    main:

    ch_versions = Channel.empty()

    //
    // Assemble reads with Unicycler
    //
    SKESA (
        reads
    )
    ch_versions = ch_versions.mix(SKESA.out.versions.first())

    //
    // Unzip scaffolds file
    //
    GUNZIP_SCAFFOLDS (
        SKESA.out.contigs
    )
    ch_versions = ch_versions.mix(GUNZIP_SCAFFOLDS.out.versions.first())

    //
    // Filter for empty scaffold files
    //
    GUNZIP_SCAFFOLDS
        .out
        .gunzip
        .filter { meta, scaffold -> scaffold.size() > 0 }
        .set { ch_scaffolds }

    //
    // Downstream assembly steps
    //
    ASSEMBLY_QC_SKESA (
        ch_scaffolds,
        fasta,
        gff,
        blast_db,
        blast_header
    )
    ch_versions = ch_versions.mix(ASSEMBLY_QC_SKESA.out.versions)

    emit:
    skesa_scaffolds          = SKESA.out.contigs                  // channel: [ val(meta), [ contigs ] ]

    blast_txt          = ASSEMBLY_QC_SKESA.out.blast_txt          // channel: [ val(meta), [ txt ] ]
    blast_filter_txt   = ASSEMBLY_QC_SKESA.out.blast_filter_txt   // channel: [ val(meta), [ txt ] ]

    quast_results      = ASSEMBLY_QC_SKESA.out.quast_results      // channel: [ val(meta), [ results ] ]
    quast_tsv          = ASSEMBLY_QC_SKESA.out.quast_tsv          // channel: [ val(meta), [ tsv ] ]

    plasmidid_html     = ASSEMBLY_QC_SKESA.out.plasmidid_html     // channel: [ val(meta), [ html ] ]
    plasmidid_tab      = ASSEMBLY_QC_SKESA.out.plasmidid_tab      // channel: [ val(meta), [ tab ] ]
    plasmidid_images   = ASSEMBLY_QC_SKESA.out.plasmidid_images   // channel: [ val(meta), [ images/ ] ]
    plasmidid_logs     = ASSEMBLY_QC_SKESA.out.plasmidid_logs     // channel: [ val(meta), [ logs/ ] ]
    plasmidid_data     = ASSEMBLY_QC_SKESA.out.plasmidid_data     // channel: [ val(meta), [ data/ ] ]
    plasmidid_database = ASSEMBLY_QC_SKESA.out.plasmidid_database // channel: [ val(meta), [ database/ ] ]
    plasmidid_fasta    = ASSEMBLY_QC_SKESA.out.plasmidid_fasta    // channel: [ val(meta), [ fasta_files/ ] ]
    plasmidid_kmer     = ASSEMBLY_QC_SKESA.out.plasmidid_kmer     // channel: [ val(meta), [ kmer/ ] ]

    versions           = ch_versions                        // channel: [ versions.yml ]
}


