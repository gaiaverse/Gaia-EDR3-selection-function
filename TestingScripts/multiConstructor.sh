#!/bin/bash
addqueue -q bigmem -c "3 hours" -n 1 -m 30 -o flat.txt ./launchConstructor.sh 0
addqueue -q bigmem -c "3 hours" -n 1 -m 30 -o gaps.txt ./launchConstructor.sh 1
addqueue -q bigmem -c "3 hours" -n 1 -m 30 -o galaxy.txt ./launchConstructor.sh 2
addqueue -q bigmem -c "3 hours" -n 1 -m 30 -o full.txt ./launchConstructor.sh 3
