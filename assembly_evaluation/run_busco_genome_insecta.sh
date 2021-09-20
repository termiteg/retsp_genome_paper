#!/bin/sh
#$ -S /bin/bash

cd ${PBS_O_WORKDIR}

MODE=genome
INPUT=References/final.assembly.fasta

LINEAGE=insecta
# eukaryota 
NCPU=24

OUTPUT=busco_out_`basename ${INPUT}`_${LINEAGE}

busco -m $MODE -i $INPUT -o $OUTPUT -l $LINEAGE -c $NCPU

