#!/bin/bash

#check arguments
if [ -z "$1" ]; then
	echo "Must input p-value limit"
	exit
fi

if [ -z "$2" ]; then
	echo "No directory inputted. Using qsub_output as analysis directory."
	output_dir=qsub_output 
else
	output_dir=$2
fi

limit=$1
output_file=${output_dir}/processed_${limit/./}.txt
if [ -f "$output_file" ]; then output_file=${output_file/.txt/}_1.txt; fi

for s in {1489..1679}; do 
	pval=`cat ${output_dir}/error_${s}.txt | sed -n 8p | sed -e 's/[eE]+*/\\*10\\^/'`
	if (( $(echo "${limit} > ${pval}" | bc -l) )); then 
		direction=`cat ${output_dir}/error_${s}.txt | sed -n 9p` 
		nam=`cat ${output_dir}/error_${s}.txt | sed -n 3p` 
		echo "${pval} : ${direction} : ${nam}" >> ${output_file}
	fi
done
