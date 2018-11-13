#!/bin/bash

export PREFIX="/home/tomswinburne/.local"

export LD_LIBRARY_PATH=$PREFIX"/lib":$LD_LIBRARY_PATH
export PATH="/home/tomswinburne/miniconda2/bin":$PREFIX"/bin":$PATH

NP=4
echo "localhost slots="$NP > hostfile # just in case...

rm -rf ./pafi

ln -s ../build/pafi pafi

mpirun --hostfile hostfile -np ${NP} ./pafi 

