#!/bin/bash
source /broad/software/scripts/useuse
reuse UGER

# $1 = job name
# $2 = script
# $3 = potential argument

#qsub -N $1 -l h_vmem=32g -l h_rt=24:00:00 -b y -p -10 -cwd -j y -V $2 $3

qsub -N $1 -pe smp 8 -binding linear:8 -l h_rt=24:00:00 -b y -p -10 -cwd -j y -V $2 $3

