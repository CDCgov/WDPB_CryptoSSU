#!/bin/bash

source /etc/profile

#$ -q short.q
#$ -N mlst-jacc
#$ -o mlst-jacc.out
#$ -e mlst-jacc.err
#$ -pe smp 1
#$ -l h_rt=2:00:00
#$ -l h_vmem=4G
#$ -cwd
#$ -V

GENOME=`realpath $1`
LIST=`realpath $2`
PREFIX=$3

# Activate Conda env for python
# module load miniconda3
# source activate /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/MEPI/ngr8/my_conda_envs/apicomplexa
module load Python3
# Run mlst_jaccard script in 1-v-all mode.
for GENOME2 in `cat $LIST`; do \
    python /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/MEPI/ngr8/bin/mlst_jaccard.py \
        -f1 $GENOME \
        -f2 $GENOME2
done > $PREFIX.mlst.tab

# Note: this script was originally written to function with Etoki's output MLST fasta files. 
# Future maintenance may be needed to make it more compatible with other formats.

# Column headers from output tab file:
# 1: Genome1
# 2: Genome2
# 3: Union of gene sets (total number of genes present in both genomes)
# 4: Intersection of gene sets (shared genes between both genomes)
# 5: Gene content difference, aka Union (minus) intersection
# 6: Jaccard *similarity* index, aka Intersection (divided by) union
# 7: Shared alleles, computed only from shared genes in Intersection (does not count genes present in only one genome as an allele difference)
# 8: Allele jaccard *similarity*, aka shared alleles (divided by) total alleles in Intersection