#!/bin/bash

# generate a list of genome FASTA files (unzipped)
# realpath *.whatever.fasta > genomes.list
# Run Etoki qsub script for each genome in the list file we just created
for i in `cat /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/Outbreak/NE_093024_latest/genomelist.txt`; do \
    PREFIX=$(basename $i .scaffolds.fa)
    qsub /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/MEPI/ngr8/bin/etoki.qsub.sh \
        $i \
        /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/MEPI/ngr8/assets/reference_alleles.fasta \
        /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/MEPI/ngr8/assets/etoki_index_md5hash.csv \
        /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/Outbreak/NE_093024_latest/Etoki_results/$PREFIX
done
# Run pairwise MLSt comparisons for each pair of Etoki outputs
# realpath *.etoki.fasta > etoki.list

for i in `cat /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/CDS_cryptoSSU_result/CDS/Outbreak_OH-IA-IL-SC-TX_8-15-24_081524/samplesheet_IIa_UPDATED/Etoki.list`; do \
    PREFIX=$(basename $i .etoki.fasta)
    qsub /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/MEPI/ngr8/bin/mlst_jaccard.qsub.sh \
        $i \
        Etoki.list \
        $PREFIX
done

# Cat all *.mlst.tab files together from outputs of pairwise MLST comparisons
# cat *.mlst.tab > All.etoki_pairs.tab

# Column headers from output `All.etoki_pairs.tab` file:
# 1: Genome1
# 2: Genome2
# 3: Union of gene sets (total number of genes present in both genomes)
# 4: Intersection of gene sets (shared genes between both genomes)
# 5: Gene content difference, aka Union (minus) intersection
# 6: Jaccard *similarity* index, aka Intersection (divided by) union
# 7: Shared alleles, computed only from shared genes in Intersection (does not count genes present in only one genome as an allele difference)
# 8: Allele jaccard *similarity*, aka shared alleles (divided by) total alleles in Intersection
