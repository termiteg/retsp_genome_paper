#!/bin/bash
set -e
set -u
set -o pipefail

BISMARK_BASE=~/bio/applications/bismark_v0.16.1
GENOME_DIR=../../../Data/TermiteG/Sequences/Rspe02
BAM=bismark_150618x6.merged.sorted.bam
NCPU=20
MEM=100G
OUTDIR=bismark_methylcall_result_150618x6_merged


ulimit -a

##
# single end reads

mkdir $OUTDIR

$BISMARK_BASE/bismark_methylation_extractor \
  -s \
  --bedGraph \
  --cytosine_report \
  --CX \
  --counts \
  --multicore $NCPU \
  --buffer_size $MEM \
  --scaffolds \
  --genome_folder $GENOME_DIR \
  --ignore 4 \
  -o $OUTDIR \
  $BAM

