#!bin/bash

for i in {3001..4000}
do
	echo ""
	echo ""
	echo "BEGINNING TEST RUN WITH RANDOM SEED "$i
	FILE="Output/Test"$i
	mpirun -n 1 ./gaia -f $FILE -s $i
done
