#!/bin/bash
FILE=$1
SUFFIX=$2
DIVISION=$3
split -l$((`wc -l < $FILE$SUFFIX`/$DIVISION)) $FILE$SUFFIX $FILE"_" --additional-suffix="$SUFFIX" -a 1
