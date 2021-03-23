#!/bin/sh
#$ -S /bin/sh
#$ -q small
#$ -pe smp 10
#$ -cwd
 
/home/yhaya/pkg/subread-1.5.0-Linux-x86_64/bin/featureCounts \
-a ../RsGM8.gtf \
-o ../fcounts/featureCounts.out \
-s 2 -T 10 \
../tophat/afb1/accepted_hits.bam \
../tophat/afb2/accepted_hits.bam \
../tophat/afb3/accepted_hits.bam \
../tophat/afh1/accepted_hits.bam \
../tophat/afh2/accepted_hits.bam \
../tophat/afh3/accepted_hits.bam \
../tophat/amb1/accepted_hits.bam \
../tophat/amb2/accepted_hits.bam \
../tophat/amb3/accepted_hits.bam \
../tophat/amh1/accepted_hits.bam \
../tophat/amh2/accepted_hits.bam \
../tophat/amh3/accepted_hits.bam \
../tophat/sfb1/accepted_hits.bam \
../tophat/sfb2/accepted_hits.bam \
../tophat/sfb3/accepted_hits.bam \
../tophat/sfh1/accepted_hits.bam \
../tophat/sfh2/accepted_hits.bam \
../tophat/sfh3/accepted_hits.bam \
../tophat/smb1/accepted_hits.bam \
../tophat/smb2/accepted_hits.bam \
../tophat/smb3/accepted_hits.bam \
../tophat/smh1/accepted_hits.bam \
../tophat/smh2/accepted_hits.bam \
../tophat/smh3/accepted_hits.bam \
../tophat/wfb1/accepted_hits.bam \
../tophat/wfb2/accepted_hits.bam \
../tophat/wfb3/accepted_hits.bam \
../tophat/wfh1/accepted_hits.bam \
../tophat/wfh2/accepted_hits.bam \
../tophat/wfh3/accepted_hits.bam \
../tophat/wmb1/accepted_hits.bam \
../tophat/wmb2/accepted_hits.bam \
../tophat/wmb3/accepted_hits.bam \
../tophat/wmh1/accepted_hits.bam \
../tophat/wmh2/accepted_hits.bam \
../tophat/wmh3/accepted_hits.bam
