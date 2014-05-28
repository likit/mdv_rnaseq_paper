#!/bin/sh -login
#PBS -l nodes=1:ppn=5,mem=24gb,walltime=24:00:00
#PBS -m abe
#PBS -M preeyano@msu.edu
#PBS -N InterproScan_${PBS_JOBID}

module load iprscan/5.44.0

cd ${PBS_O_WORKDIR} 

iprscan -i ${input}
