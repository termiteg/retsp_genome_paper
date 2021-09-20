INPUT={{FASTQ_FILE}}
 
QV=30
QB=33
 
ADPT1=GTCATAGCTGTTTCCTGCTGGCCGTCGTTTTAC
 
INPUT_EXT=`basename $INPUT .gz`
OUTPUT=`basename $INPUT_EXT .fastq`.clnq$QV.fq
 
cutadapt -a $ADPT1 --overlap=8 --quality-base=$QB -q $QV $INPUT  > $OUTPUT
