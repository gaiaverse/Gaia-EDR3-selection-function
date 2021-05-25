#!bin/bash

mpirun -n 20 ./gaia -f Output/Diagnostic3_mum0_mut0_sigmat3_lt36_space54/ -t ../../Data/ShuffledData/ -s 3 -l 500 -g 1e-5 > Output/Testing.txt
mpirun -n 20 ./gaia -f Output/Diagnostic4_mum0_mut0_sigmat3_lt36_space65/ -t ../../Data/ShuffledData/ -s 3 -l 500 -g 1e-5 > Output/Testing.txt
mpirun -n 20 ./gaia -f Output/Diagnostic5_mum0_mut0_sigmat3_lt36_space76/ -t ../../Data/ShuffledData/ -s 3 -l 500 -g 1e-5 > Output/Testing.txt
mpirun -n 20 ./gaia -f Output/Diagnostic6_mum1_mut0_sigmat3_lt36_space65/ -t ../../Data/ShuffledData/ -s 3 -l 500 -g 1e-5 > Output/Testing.txt
mpirun -n 20 ./gaia -f Output/Diagnostic7_mum2_mut0_sigmat3_lt36_space65/ -t ../../Data/ShuffledData/ -s 3 -l 500 -g 1e-5 > Output/Testing.txt
mpirun -n 20 ./gaia -f Output/Diagnostic8_mum25_mut0_sigmat3_lt36_space65/ -t ../../Data/ShuffledData/ -s 3 -l 500 -g 1e-5 > Output/Testing.txt
