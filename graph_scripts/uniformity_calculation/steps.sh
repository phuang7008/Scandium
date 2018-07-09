#!/bin/bash

# this script is used to setup steps for the calculation of the uniformity scores
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/steps.sh Peak_size version"
	exit
fi

path=`pwd`
peak_size=$1
version=$2

echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/uniformity.sh $path $peak_size $version" | msub -V -d $path -q high -A proj-dm0001 -j oe -N uniformity -l nodes=1:ppn=1,mem=5gb
