#!/bin/bash

# this script is used to setup steps for the calculation of the uniformity scores
#
if [[ $# -ne 3 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/projects/group_related_works/Novaseq_bias/uniformity.sh file_path Peak_size version"
	exit
fi

path=$1
peak_size=$2
version=$3

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

	echo "$bambasename"

	echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison.pl $bFile $version $peak_size >> uniformity_scores"

	/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison.pl $bFile $version $peak_size >> uniformity_scores

	#break
done
