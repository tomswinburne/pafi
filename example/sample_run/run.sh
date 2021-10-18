# create / replace soft link
touch pafi
rm pafi
ln -s ../../build/pafi pafi

# create dump folder
mkdir -p dumps

NP=1

# hostfile not always required
echo "localhost slots="${NP} > hostfile

mpirun --hostfile hostfile -np ${NP} ./pafi

# clean up
rm hostfile
