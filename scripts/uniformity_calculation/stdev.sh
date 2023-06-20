#!/bin/bash

# this script is used to setup steps for the calculation of the uniformity scores
#
if [[ $# -ne 1 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/projects/group_related_works/Novaseq_bias/stdev_steps.sh max_outlier"
	exit
fi

max_outlier=$1

path=`pwd`
echo "$path"

id=1

for bFile in `ls $path/*WGS_between1x_150x_REPORT.txt`
#for bFile in `ls $path/*wo_centromeres`
do
	bambasename=""
	if [[ $bFile = *"_wo_centromeres"* ]]; then
		bambasename=$(basename "$bFile" _wo_centromeres) 
	fi

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

	#echo "$bambasename"

	echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/sequence_stats_calculation.pl $bFile $max_outlier " | msub -V -d $path -q high -A proj-dm0001 -j oe -N stdev_$id -l nodes=1:ppn=1,mem=5gb

	id=$((id+1))

done
