!/bin/sh
#$ -S /bin/sh
#$ -q small
#$ -cwd
 
cd /home/yhaya/2014_Rspe_Genome/151130_InParanoid/multiparanoid
 
/home/yhaya/pkg/multiparanoid/multiparanoid.pl -unique 1 \
-species Znev+Mnat+Rspe
 
