#!/bin/sh

if [ -z "$1" ]; then
	echo "Must input p-value limit"
	exit
fi

if [ "$2" != "" ]; then
	output_file=$2
	echo "Writing to ${output_file}"
fi

filename=drug_coxph_summ.txt

if [ -z ${filename} ]; then
	echo "Must have ${filename} in this directory"
	exit
fi

limit=$1

if [ -z "${output_file}" ]; then
	echo "Checking for limit of ${limit}\n"
else
	echo "Checking for limit of ${limit}" >> ${output_file}
fi
grep -wn "DOXAZOSIN_dur" ${filename} | while read -r line ; do 
	var=`echo ${line} | awk -F " " '{print $NF}' | sed -e 's/[eE]+*/\\*10\\^/'`
	num=`echo ${line} | awk -F ":" '{print $1}'`
	dir=`echo ${line} | awk -F " " '{print $2}'`
	if [ "${var}" != "+" ]; then #doesn't check if numeric
		if (( $(echo "${limit} > ${var}" | bc -l) )); then 
			wanted=$(( $num-7 ))
			if [ -z "${output_file}" ]; then
				echo "${var}:${dir}:`cat ${filename} | sed -n ${wanted}p`"
			else
				echo "${var}:${dir}:`cat ${filename} | sed -n ${wanted}p`" >> ${output_file}
			fi
		fi
	fi
done
