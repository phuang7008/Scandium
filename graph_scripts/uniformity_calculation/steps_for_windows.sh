#!/bin/bash

# this script is used to setup steps for the calculation of the uniformity scores
#
if [[ $# -ne 1 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/steps_for_windows.sh version"
	exit
fi

path=`pwd`
version=$1

for bFile in `ls $path/*WGS_between1x_150x_REPORT.txt`
#for bFile in `ls $path/*bam*WGS_between1x_150x_REPORT.txt`
#for bFile in `ls $path/*cram*WGS_between1x_150x_REPORT.txt`
#for bFile in `ls $path/*hgv.WGS_between1x_150x_REPORT.txt`
#for bFile in `ls $path/*wo_centromeres`
do
	echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/window_size_uniformity.sh $bFile $version" | msub -V -d $path -q high -A proj-dm0001 -j oe -N uniformity_windows -l nodes=1:ppn=1,mem=5gb

done
