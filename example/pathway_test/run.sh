# create / replace soft link
touch pafi
rm pafi
ln -s ../../build/pafi pafi

NP=1

# create dump folder
mkdir -p dumps

# hostfile not always required
echo "localhost slots="${NP} > hostfile

mpirun --hostfile hostfile -np ${NP} ./pafi

# clean up
rm hostfile
