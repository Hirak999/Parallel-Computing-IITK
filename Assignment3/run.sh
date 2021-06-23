#!/bin/bash

make

~/UGP/allocator/src/allocator.out 4 4
mpiexec -np 1  -hostfile hosts ./src.x tdata.csv
mpiexec -np 2  -hostfile hosts ./src.x tdata.csv
mpiexec -np 4  -hostfile hosts ./src.x tdata.csv
~/UGP/allocator/src/allocator.out 2 1
mpiexec -np 2 -hostfile hosts ./src.x tdata.csv
~/UGP/allocator/src/allocator.out 4 2
mpiexec -np 4 -hostfile hosts ./src.x tdata.csv
~/UGP/allocator/src/allocator.out 8 4
mpiexec -np 8 -hostfile hosts ./src.x tdata.csv

