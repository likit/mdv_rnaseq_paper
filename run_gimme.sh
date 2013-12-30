#!/bin/sh -login
#PBS -l nodes=1:ppn=1,mem=32gb,walltime=24:00:00
#PBS -M preeyano@msu.edu
#PBS -m abe
#PBS -N Gimme_${PBS_JOBID}

module load bxPython
module load pygr
module load matplotlib

cd /mnt/ls12/preeyanon/gal4models_new
python ~/gimme/src/gimme.py -r ${ref} ${input} > ${input}.bed 2>${input}.gimme.log
