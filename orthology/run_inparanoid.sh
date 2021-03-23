#!/bin/sh
#$ -S /bin/sh
#$ -q medium
#$ -pe smp 2
#$ -t 1-435
#$ -tc 80
#$ -cwd
 
export PATH=/home/yhaya/pkg/blast-2.2.26/bin:$PATH
ulimit -s unlimited
 
array1=(`ls -1 ../seqfiles`)
array2=(${array1[@]})
 
pairs=()
for i in ${array1[@]}
do
    array2=("${array2[@]:1}")
    for j in  ${array2[@]}
    do
        p="${i}-${j}"
        pairs+=( $p )
    done
done
 
index=`expr ${SGE_TASK_ID} - 1`
pair=${pairs[index]}
SP1=${pair%-*}
SP2=${pair#*-}
 
mkdir ../inparanoid/$pair
cd ../inparanoid/$pair
 
cp /home/yhaya/pkg/inparanoid_4.1/*.pl \
  /home/yhaya/pkg/inparanoid_4.1/seqstat.jar \
  /home/yhaya/pkg/inparanoid_4.1/BLOSUM62 \
  /home/yhaya/pkg/inparanoid_4.1/PAM* \
  ../../seqfiles/$SP1 \
  ../../seqfiles/$SP2 \
  ./
 
./inparanoid_2threads.pl \
$SP1 $SP2
 
cp sqltable.* ../../multiparanoid
rm $SP1 $SP2 PAM* BLOSUM62 seqstat.jar *.pl
