#! /bin/bash 

InputSNP=$1
InputTF=$2
Threshold=$3


while read -r line
do
	count=$(grep -i "${line}" ${InputSNP} | wc -l)
	if [ "$count" -ne 0 ]; then
		echo "${line} being Run"
		#echo $count
		#half_count=$((count / 2))
		#echo $half_count
		mkdir -p "TF/${line}"
		grep -i "${line}" ${InputSNP} > "TF/${line}/SNP.tsv"
		#grep -i "${line}" ../mutations/FalseABS.txt | tail -"${half_count}" >> "TF/${line}/SNP.tsv"
		python execute.py -t "$line" -p $Threshold
	fi
done < $InputTF

