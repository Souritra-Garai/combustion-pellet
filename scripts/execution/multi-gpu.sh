#!/bin/bash

usage="Usage: $0 executable inp1 inp2 inp3 ... inpN"

if [ "$#" == "0" ]; then
	echo "$usage"
	exit 1
fi

assigned_gpu=0
gpus_count=$(nvidia-smi -L | wc -l)
executable="$1"
shift

while (($#)); do
	job=$1
	echo "Run ${job} on ${assigned_gpu}"
	export CUDA_VISIBLE_DEVICES=${assigned_gpu}
	export WSLENV=CUDA_VISIBLE_DEVICES/w
	
	${executable} --phi ${job} &

	assigned_gpu=$(($assigned_gpu + 1))
	
	if [ $assigned_gpu -eq $gpus_count ]; then
		wait
	fi

	assigned_gpu=$(($assigned_gpu % $gpus_count))
	shift
done

wait