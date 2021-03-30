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

for file in `ls $path/*WGS_uniformity_REPORT.txt`
do
    filename=$(basename -- "$file")
    filename="${filename%.*}"
    filename="${filename%.*}"
    filename="${filename%.*}"
    filename="${filename%.*}"
	echo "$filename"

	echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison.pl $file $version $peak_size >> uniformity_scores"

	/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison.pl $file $version $peak_size >> uniformity_scores

	#break
done
