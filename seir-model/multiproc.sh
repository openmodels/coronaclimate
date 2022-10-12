#!/bin/bash

if [ "$#" -ne 1 ]; then
    for i in $(seq 1 $2); do
	nohup Rscript $1 "{@:3}" > /dev/null 2>&1 &
	sleep 5
    done
else
    Rscript $1 "{@:3}"
fi
