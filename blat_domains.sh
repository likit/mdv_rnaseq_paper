#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=24gb,walltime=4:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Blat_domains_${PBS_JOBID}

cd ${PBS_O_WORKDIR} 

blat ${index} ${input} -q=prot -t=dnax -noHead ${input}.psl
