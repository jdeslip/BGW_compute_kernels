#!/bin/bash

#export KMP_AFFINITY=scatter
export KMP_AFFINITY=compact
get_hostfile

#nthreads="1 2 4 8 16 32"
nthreads="32 16 8 4 2"
for nthread in $nthreads; do
	echo $nthread
	export OMP_NUM_THREADS=$nthread
	mpirun -n 1 -hostfile hostfile.$PBS_JOBID -ppn 1 ${PWD}/kernel_BSE.x 2048  6 15 2  2 8 20 &> host_${nthread}.log
done

