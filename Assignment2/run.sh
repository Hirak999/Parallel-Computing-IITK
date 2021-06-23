#!/bin/bash
#!/usr/bin/env python3.1

chmod +x script.py
python3 ./script.py 4 4 8

make

for (( j=1; j<=10; j++ ))
do
{
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 2048 1 
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 32768 1 
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 262144 1 
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 2048 1
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 32768 1 
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 262144 1 
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 2048 1
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 32768 1 
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 262144 1 
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 2048 1 
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 32768 1 
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 262144 1 
}
done

for (( j=1; j<=10; j++ ))
do
{
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 2048 2
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 32768 2
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 262144 2
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 2048 2
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 32768 2
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 262144 2
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 2048 2
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 32768 2
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 262144 2
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 2048 2
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 32768 2
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 262144 2
}
done

for (( j=1; j<=10; j++ ))
do
{
mpiexec -np 4 -ppn 1 -f hostfile ./src.x  2048 3
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 32768 3
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 262144 3
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 2048 3
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 32768 3
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 262144 3
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 2048 3
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 32768 3
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 262144 3
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 2048 3
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 32768 3
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 262144 3
}
done

for (( j=1; j<=10; j++ ))
do
{
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 2048 4
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 32768 4
mpiexec -np 4 -ppn 1 -f hostfile ./src.x 262144 4
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 2048 4
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 32768 4
mpiexec -np 32 -ppn 8 -f hostfile ./src.x 262144 4
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 2048 4
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 32768 4
mpiexec -np 16 -ppn 1 -f hostfile ./src.x 262144 4
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 2048 4
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 32768 4
mpiexec -np 128 -ppn 8 -f hostfile ./src.x 262144 4
}
done

python3 plot.py
