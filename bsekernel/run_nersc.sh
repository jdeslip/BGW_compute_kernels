#!/bin/bash
export OMP_NUM_THREADS=240
#export OMP_NUM_THREADS=236
#export OMP_NUM_THREADS=118
#export OMP_NUM_THREADS=5
export KMP_AFFINITY=balanced,verbose
#export KMP_PLACE_THREADS=59C,4T
get_micfile
mpirun.mic -n 1 -hostfile micfile.$PBS_JOBID -ppn 1 ${PWD}/kernel_BSE.x 2048  6 15 2  2 8 20
#multiply time by 8!!
#mpirun.mic -n 1 -hostfile micfile.$PBS_JOBID -ppn 1 ${PWD}/kernel_BSE.x 2048  6 15 2  2 8 10
