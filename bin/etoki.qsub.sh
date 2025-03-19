#!/bin/bash

source /etc/profile

#$ -q short.q
#$ -N etoki
#$ -o etoki.out
#$ -e etoki.err
#$ -pe smp 1
#$ -l h_rt=3:00:00
#$ -l h_vmem=4G
#$ -cwd
#$ -V

FASTA=`realpath $1`
REF=`realpath $2`
DB=`realpath $3`
PREFIX=$4

# Run Etoki
singularity exec /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/MEPI/ngr8/dev/wdpb_vapor1_mlstknockout/bin/etoki_latest.sif \
    EToKi.py MLSType \
        -i $FASTA \
        -r $REF \
        -k $PREFIX \
        -d $DB \
        -o $PREFIX.etoki.fasta

# Default values for Crypto
# REF == /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/MEPI/ngr8/assets/reference_alleles.fasta
# DB  == /scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/MEPI/ngr8/assets/etoki_index_md5hash.csv
