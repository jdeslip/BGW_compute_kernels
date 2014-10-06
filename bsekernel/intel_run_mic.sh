#!/bin/bash
export OMP_NUM_THREADS=240
export KMP_AFFINITY=compact,verbose
/home/mic/kernel_BSE.x 2048  6 15 2  2 8 20
