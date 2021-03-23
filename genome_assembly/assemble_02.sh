#!/bin/sh

#$ -S /bin/bash
#$ -l nc=24
#$ -cwd

export PATH=/home/shige/bio/Applications/allpathslg-47878/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/gcc-4.7.3/lib64::$LD_LIBRARY_PATH


ulimit -s 100000

RunAllPathsLG \
 PRE=$PWD \
 REFERENCE_NAME=Rspe02 \
 DATA_SUBDIR=data \
 RUN=run \
 OVERWRITE=True \
 VAPI_WARN_ONLY=True \
 THREADS=30 \
 > assemble_02.out

