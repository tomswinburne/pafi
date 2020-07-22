export PREFIX=$HOME"/.local"
export LD_LIBRARY_PATH=$PREFIX"/lib":$LD_LIBRARY_PATH
export PATH=$PREFIX"/bin":$PATH




cat preconfig.xml | sed "s/CC/$1/g" > config.xml

mpirun -np $1 ./pafi

#touch pafi
#rm pafi
#ln -s ../build/pafi pafi

#mkdir -p dumps

#NP=4
## test example, not always required
#echo "localhost slots="$NP > hostfile

#mpirun --hostfile hostfile -np ${NP} ./pafi
