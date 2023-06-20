#!/bin/bash

# this script is used to setup steps for the calculation of the uniformity scores
#
if [[ $# -ne 1 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/projects/group_related_works/Novaseq_bias/subtract_centromeres.sh centromere_file"
	exit
fi

centromeres=$1

path=`pwd`
resultPath=$path/no_centromeres
echo "$path"

for bFile in `ls $path/*WGS_between1x_150x_REPORT.txt`
do
	bambasename=""
	if [[ $bFile = *"hgv.cram"* ]]; then
		bambasename=$(basename "$bFile" .hgv.cram.WGS_between1x_150x_REPORT.txt) 
	fi

	if [[ $bFile = *"hgv.bam"* ]]; then 
		bambasename=$(basename "$bFile" .hgv.bam.WGS_between1x_150x_REPORT.txt) 
	fi

	if [[ $bFile = *"realigned.recal.bam"* ]]; then
		bambasename=$(basename "$bFile" .realigned.recal.bam.WGS_between1x_150x_REPORT.txt) 
	fi

	if [[ $bFile = *"realigned.recal.cram"* ]]; then
		bambasename=$(basename "$bFile" .realigned.recal.cram.WGS_between1x_150x_REPORT.txt) 
	fi

	if [[ $bFile = *"cram.WGS"* ]]; then
		bambasename=$(basename "$bFile" .cram.WGS_between1x_150x_REPORT.txt)
	fi

	#echo "$bambasename"

	echo "less $bFile | grep -v '#' | /hgsc_software/BEDTools/latest/bin/bedtools subtract -a - -b $centromeres > $resultPath/${bambasename}_wo_centromeres" | msub -V -d $resultPath -q high -A proj-dm0001 -j oe -N centromere -l nodes=1:ppn=1,mem=2gb

	#break
done
