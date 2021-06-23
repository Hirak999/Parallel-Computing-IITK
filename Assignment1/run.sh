#!/bin/bash


~/UGP/allocator/src/allocator.out 64 8

make

mpiexec -np 16 -hostfile hosts ./src.x 256 50
mpiexec -np 16 -hostfile hosts ./src.x 1024 50
mpiexec -np 16 -hostfile hosts ./src.x 4096 50
mpiexec -np 16 -hostfile hosts ./src.x 16384 50
mpiexec -np 16 -hostfile hosts ./src.x 65536 50
mpiexec -np 16 -hostfile hosts ./src.x 262144 50
mpiexec -np 16 -hostfile hosts ./src.x 1048576 50
mpiexec -np 36 -hostfile hosts ./src.x 256 50
mpiexec -np 36 -hostfile hosts ./src.x 1024 50
mpiexec -np 36 -hostfile hosts ./src.x 4096 50
mpiexec -np 36 -hostfile hosts ./src.x 16384 50
mpiexec -np 36 -hostfile hosts ./src.x 65536 50
mpiexec -np 36 -hostfile hosts ./src.x 262144 50
mpiexec -np 36 -hostfile hosts ./src.x 1048576 50
mpiexec -np 49 -hostfile hosts ./src.x 256 50
mpiexec -np 49 -hostfile hosts ./src.x 1024 50
mpiexec -np 49 -hostfile hosts ./src.x 4096 50
mpiexec -np 49 -hostfile hosts ./src.x 16384 50
mpiexec -np 49 -hostfile hosts ./src.x 65536 50
mpiexec -np 49 -hostfile hosts ./src.x 262144 50
mpiexec -np 49 -hostfile hosts ./src.x 1048576 50
mpiexec -np 64 -hostfile hosts ./src.x 256 50
mpiexec -np 64 -hostfile hosts ./src.x 1024 50
mpiexec -np 64 -hostfile hosts ./src.x 4096 50
mpiexec -np 64 -hostfile hosts ./src.x 16384 50
mpiexec -np 64 -hostfile hosts ./src.x 65536 50
mpiexec -np 64 -hostfile hosts ./src.x 262144 50
mpiexec -np 64 -hostfile hosts ./src.x 1048576 50

python3 plot16.py
python3 plot36.py
python3 plot49.py
python3 plot64.py
