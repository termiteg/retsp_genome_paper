#!/bin/sh
#$ -S /bin/sh
#$ -q small
#$ -t 1-36
#$ -pe smp 10
#$ -cwd
 
export PATH=/home/yhaya/pkg/bowtie2-2.2.3:/home/yhaya/pkg/samtools-1.1:$PATH
idir=/home/yhaya/Analysis/2015_Rspe_RNAseq/151130_Trimming
odir=/home/yhaya/Analysis/2015_Rspe_RNAseq/151201_TopHat_RsGM8
 
x=`ls ${idir}/*.fq.gz | xargs -I % basename % | sed -e 's/\..\+//' | uniq | head -n ${SGE_TASK_ID} | tail -n 1`
 
/home/yhaya/pkg/tophat-2.1.0.Linux_x86_64/tophat \
-p 10 -o ${odir}/tophat/${x} \
--library-type fr-firststrand \
--no-novel-juncs \
--GTF /home/yhaya/2014_Rspe_Genome/151008_GeneModelCuration/RsGM8.1/RsGM8.1.gff3 \
/home/yhaya/Analysis/2015_Rspe_RNAseq/151201_TopHat_RsGM8/index/Rspe2 \
