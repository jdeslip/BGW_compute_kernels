#!/bin/bash

export I_MPI_PIN_DOMAIN=omp
export OMP_NUM_THREADS=240
#export OMP_NUM_THREADS=236
#export OMP_NUM_THREADS=118
#export OMP_NUM_THREADS=5
export KMP_AFFINITY=balanced
#export KMP_PLACE_THREADS=59C,4T
get_micfile
nthreads="240 120 60"
for nthread in $nthreads; do
	echo $nthread
	export OMP_NUM_THREADS=$nthread
	mpirun.mic -n 1 -hostfile micfile.$PBS_JOBID -ppn 1 ${PWD}/../kernel_BSE.x 2048  6 15 2  2 8 20 &> host_${nthread}.log
done
