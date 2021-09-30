#!/bin/bash
addqueue -q bigmem -c "3 hours" -n 1 -m 30 ./launchConstructor.sh 0
addqueue -q bigmem -c "3 hours" -n 1 -m 30 ./launchConstructor.sh 1
addqueue -q bigmem -c "3 hours" -n 1 -m 30 ./launchConstructor.sh 2
addqueue -q bigmem -c "3 hours" -n 1 -m 30 ./launchConstructor.sh 3
