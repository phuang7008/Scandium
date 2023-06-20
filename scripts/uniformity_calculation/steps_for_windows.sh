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

for file in `ls $path/*WGS_uniformity_REPORT.txt`
do
	echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/repo/scandium/graph_scripts/uniformity_calculation/window_size_uniformity.sh $file $version" | msub -V -d $path -q high -A proj-dm0001 -j oe -N uniformity_windows -l nodes=1:ppn=1,mem=5gb
	#echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/repo/scandium/graph_scripts/uniformity_calculation/window_size_uniformity.sh $file hg38" | msub -V -d $path -q high -A proj-dm0001 -j oe -N uniformity_windows -l nodes=1:ppn=1,mem=5gb

done

#for file in `ls $path/*bam*bam.WGS_uniformity_REPORT.txt`
#do
#	echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/repo/scandium/graph_scripts/uniformity_calculation/window_size_uniformity.sh $file hg37" | msub -V -d $path -q high -A proj-dm0001 -j oe -N uniformity_windows -l nodes=1:ppn=1,mem=5gb

#done
