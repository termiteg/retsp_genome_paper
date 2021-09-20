#!/bin/bash

#===
# Template information:
# This is a template to run bismark for single end Illumina reads.
# Work with bismark version > v0.16.1
# date: May 13, 2016
# author: Shuji Shigenobu <shige@nibb.ac.jp>
#===

set -e
set -u
set -o pipefail


BISMARK_BASE=~/bio/applications/bismark_v0.16.1
GENOME_DIR=../../../Data/TermiteG/Sequences/Rspe02
NAME=%NAME%
SEQ1=%SEQ1%
NCPU=2
OUT=bismark_out_$NAME


$BISMARK_BASE/bismark \
  --bowtie2 \
  -p $NCPU \
  --pbat \
  --unmapped \
  -o $OUT \
  $GENOME_DIR \
  $SEQ1 \

echo "======="
echo "`date` -- completed on host: `hostname`"

