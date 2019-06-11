#!/bin/bash

# this script is used to setup steps for the calculation of the uniformity scores
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/projects/group_related_works/Novaseq_bias/window_size_uniformity.sh file_name version"
	exit
fi

path=`pwd`
bFile=$1
version=$2
echo "version is $version"

#w_size=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40)
#w_size=(1 2 3 4 5 6 7 8 9 10)

bambasename=""
if [[ $bFile = *"_wo_centromeres"* ]]; then
	bambasename=$(basename "$bFile" _wo_centromeres) 
fi

if [[ $bFile = *".hgv.cram.WGS_uniformity_REPORT.txt" ]]; then
	bambasename=$(basename "$bFile" .hgv.cram.WGS_uniformity_REPORT.txt) 
fi

if [[ $bFile = *".hgv.bam.WGS_uniformity_REPORT.txt" ]]; then 
	bambasename=$(basename "$bFile" .hgv.bam.WGS_uniformity_REPORT.txt) 
fi

if [[ $bFile = *".realigned.recal.bam.WGS_uniformity_REPORT.txt" ]]; then
	bambasename=$(basename "$bFile" .realigned.recal.bam.WGS_uniformity_REPORT.txt) 
fi

if [[ $bFile = *".realigned.recal.cram.WGS_uniformity_REPORT.txt" ]]; then
	bambasename=$(basename "$bFile" .realigned.recal.cram.WGS_uniformity_REPORT.txt) 
fi

if [[ $bFile = *".WGS_uniformity_REPORT.txt" ]]; then                                                         
	bambasename=$(basename "$bFile" .WGS_uniformity_REPORT.txt)                                               
fi

echo "$bambasename"

#for i in "${w_size[@]}"
for i in {1..50}
#for i in {7..7}
do
	#echo
	#echo $i
	#echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison.pl $bFile $version $i >> ${bambasename}_uniformity_windows.txt"

	#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison.pl $bFile $version $i >> ${bambasename}_uniformity_windows.txt

	#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison_autosome_X_Y.pl $bFile $version $i 1 >> ${bambasename}_uniformity_windows_all.txt

	#echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison_autosome_X_Y.pl $bFile $version $i 2 >> ${bambasename}_uniformity_windows_autosome_only.txt"

	/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison_autosome_X_Y.pl $bFile $version $i 2 >> ${bambasename}_uniformity_windows_autosome_only.txt

	#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison.pl $bFile $version $i >> ${bambasename}_uniformity_windows.txt

	#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/calculate_uniformity_comparison_w_kurtosis.pl $bFile $version $i >> ${bambasename}_uniformity_windows.txt
done
