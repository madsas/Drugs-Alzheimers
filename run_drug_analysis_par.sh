#!/bin/sh

#check argument
if [ -z "$1" ]; then
	echo "No dirname inputted. Using default directory of qsub_output"
	dirname=qsub_output
else
	dirname=$1
fi

#define the rows to read
startrow=1489
endrow=1679

#get dummy script name
ds=tempscript.sh

#create output directory for qsub
mkdir ${dirname}
wd="/afs/.ir.stanford.edu/users/s/a/sasidhar/naccdataset/greicius_project"
#define working directory


for i in  $(seq $startrow $endrow); do 
#for i in 1490; do 
	rm $ds
	echo "#!/bin/sh" >> $ds
	echo "R --vanilla --no-save --args ${i} < ${wd}/drug_analysis_edited_par.R" >> ${ds}
	qsub -o ${wd}/${dirname}/output_${i}.txt -e ${wd}/${dirname}/error_${i}.txt $ds
done
