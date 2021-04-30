#!/bin/bash

OUTDIR = $1

for i in *.csv;
do
	cat $i | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .01) print $0}' > $1$i
done
