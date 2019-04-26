#!/bin/bash

# Expand the tarball in this folder
#	 tar -xvf Fe_dumbell_example_data.tar.gz

# Adjust the environment variables in this file

export PREFIX=$HOME"/.local"
export LD_LIBRARY_PATH=$PREFIX"/lib":$LD_LIBRARY_PATH
export PATH=$PREFIX"/bin":$PATH

ln -s ../new_build/pafi pafi

NP=4
 # test example, not always required
echo "localhost slots="$NP > hostfile

mpirun --hostfile hostfile -np ${NP} ./pafi
