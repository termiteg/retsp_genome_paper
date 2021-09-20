#===
# Whole genome alignment using lastz/chain/net
# 01_prep_lastz.sh
#
# This script setup imput sequences in .nib format
#===

#=== configuration ===
# repeats should be masked (using such as RepeatMasker and WindowMasker)
TARGET_FASTA=Rspe02.assembly.mask.fa
QUERY_FASTA=Mnat_assembly_v1.0.fa.masked.lt1k

TARGET_DIR=Rspe
QUERY_DIR=Mnat
#===

TARGET=$TARGET_DIR
QUERY=$QUERY_DIR

## prep

mkdir $TARGET_DIR $QUERY_DIR
faSplit byName $TARGET_FASTA $TARGET_DIR/
faSplit byName $QUERY_FASTA $QUERY_DIR/

for f in $TARGET_DIR/*.fa; do faToNib -softMask $f `echo $f |sed -e s/.fa/.nib/`; done
for f in $QUERY_DIR/*.fa;  do faToNib -softMask $f `echo $f |sed -e s/.fa/.nib/`; done

fast.rb ids $TARGET_FASTA >$TARGET.list
fast.rb ids $QUERY_FASTA >$QUERY.list

faSize $TARGET_FASTA -detailed > $TARGET.sizes
faSize $QUERY_FASTA -detailed > $QUERY.sizes



