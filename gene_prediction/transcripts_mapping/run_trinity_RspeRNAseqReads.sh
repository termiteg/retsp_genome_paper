#!/bin/sh
#$ -S /bin/bash
#$ -l nc=24
#$ -cwd
 
#=== config ===
## REQUIRED
AHOME={{AHOME}}
SEQ1=$AHOME/left_all.fq
SEQ2=$AHOME/right_all.fq
OUT=./trinity_RspeHay_140217

  
## OPTIONAL PARAMETERS
MIN_KMER_COV=2
  
## SERVER CONF
NCPU=30
JM=256G
#===

TRINITY_BASE=~/bio/Applications/trinityrnaseq_r2013_08_14
  
echo original path: $PATH
ulimit -s unlimited
  
echo running on `hostname`
echo started at `date`
  
$TRINITY_BASE/Trinity.pl \
     --JM $JM \
     --seqType fq \
     --min_kmer_cov $MIN_KMER_COV \
     --left $SEQ1  --right $SEQ2 \
     --CPU $NCPU \
     --output $OUT \

