#!bin/bash

for i in {1..20}
do
	echo ""
	echo ""
	echo "BEGINNING TESTS RUN WITH BURNIN "$i
	
	for j in {3770..3800}
	do
		FILE="Output/Test_seed"$j"_burn"$i
		mpirun -n 1 ./gaia -f $FILE -s $j -b $i
	done
done
