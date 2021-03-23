GENOME=Rspe02.final.assembly.fasta
GFF=EVM8.evm.gff3
#EVM7.evm.out.gff3
OUTBASE=EVM7b


/home/shige/bio/Applications/EVM_r2012-06-25/EvmUtils/gff3_file_to_proteins.pl \
  $GFF    \
  $GENOME \
  prot    \
  > $OUTBASE.pep


/home/shige/bio/Applications/EVM_r2012-06-25/EvmUtils/gff3_file_to_proteins.pl \
  $GFF    \
  $GENOME \
  CDS     \
  > $OUTBASE.cds
