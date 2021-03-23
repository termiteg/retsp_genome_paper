#!/bin/sh

export PATH=/home/shige/bio/Applications/allpathslg-47878/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/gcc-4.7.3/lib64::$LD_LIBRARY_PATH

# ALLPATHS-LG needs 100 MB of stack space.  In 'csh' run 'limit stacksize 100000'.
ulimit -s 100000

mkdir -p Rspe02/data

# NOTE: The option GENOME_SIZE is OPTIONAL. 
#       It is useful when combined with FRAG_COVERAGE and JUMP_COVERAGE 
#       to downsample data sets.
#       By itself it enables the computation of coverage in the data sets 
#       reported in the last table at the end of the preparation step. 

# NOTE: If your data is in BAM format you must specify the path to your 
#       picard tools bin directory with the option: 
#
#       PICARD_TOOLS_DIR=/your/picard/tools/bin

PrepareAllPathsInputs.pl \
 DATA_DIR=$PWD/Rspe02/data \
 PLOIDY=2 \
 IN_GROUPS_CSV=in_groups2.csv \
 IN_LIBS_CSV=in_libs2.csv \
 GENOME_SIZE=1000000000 \
 OVERWRITE=False \
 VERBOSE=2 \
 > prepare_02.out 


