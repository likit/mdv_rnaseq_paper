#!/bin/sh -login
#PBS -l nodes=1:ppn=4,mem=24gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N CDHIT_${PBS_JOBID}

export PATH=$PATH:/opt/software/CDHIT/4.5.6--GCC-4.4.5

cd ${PBS_O_WORKDIR}
cd-hit-est -T 4 -d 0 -c 1.0 -M 24000 -i all.fa.clean -o all.fa.clean.nr
