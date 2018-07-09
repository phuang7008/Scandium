#!/bin/bash

path=`pwd`
id=1
size=1000
version=$1

## declare an array variable
#declare -a chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")
declare -a chromosomes=("12" "15" "19" "Y")

for chrom in "${chromosomes[@]}"
do
	for file in `ls $path/*between1x_150x*`
	do
		echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_high_cov_graph.R $file $chrom 1 $size $version " | msub -q high -l nodes=1:ppn=1,mem=32gb -V -d $path -N graph$id -o out$id -e err$id -A proj-dm0001

		id=$((id+1))
	done
done
