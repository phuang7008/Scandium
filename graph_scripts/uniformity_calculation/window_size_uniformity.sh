#!/bin/bash

# this script is used to setup steps for the calculation of the uniformity scores
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/projects/group_related_works/Novaseq_bias/window_size_uniformity.sh file_name version"
	exit
fi

path=`pwd`
file=$1
version=$2
echo "version is $version"

filename=$(basename -- "$file")
filename="${filename%.*}"
filename="${filename%.*}"
filename="${filename%.*}"
filename="${filename%.*}"

#for i in "${w_size[@]}"
for i in {1..50}
#for i in {7..7}
do
	#echo
	#echo $i
	#echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/repo/scandium/graph_scripts/uniformity_calculation/calculate_uniformity_comparison.pl $file $version $i >> ${filename}_uniformity_windows.txt"

	#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/repo/scandium/graph_scripts/uniformity_calculation/calculate_uniformity_comparison.pl $file $version $i >> ${filename}_uniformity_windows.txt

	#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/repo/scandium/graph_scripts/uniformity_calculation/calculate_uniformity_comparison_autosome_X_Y.pl $file $version $i 1 >> ${filename}_uniformity_windows_all.txt

	#echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/repo/scandium/graph_scripts/uniformity_calculation/calculate_uniformity_comparison_autosome_X_Y.pl $file $version $i 2 >> ${filename}_uniformity_windows_autosome_only.txt"

	/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/repo/scandium/graph_scripts/uniformity_calculation/calculate_uniformity_comparison_autosome_X_Y.pl $file $version $i 2 >> ${filename}_uniformity_windows_autosome_only.txt

	#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/repo/scandium/graph_scripts/uniformity_calculation/calculate_uniformity_comparison.pl $file $version $i >> ${filename}_uniformity_windows.txt

	#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/repo/scandium/graph_scripts/uniformity_calculation/calculate_uniformity_comparison_w_kurtosis.pl $file $version $i >> ${filename}_uniformity_windows.txt
done
