#!/bin/bash

path=`pwd`
id=1
size=200
chrom_id=7

#for file in `ls $path/*WGS*between*_REPORT.txt`
#for file in `ls $path/*cram*between1x_1000x_REPORT.txt`
for file in `ls $path/*wo_centromeres`
do 

	# for hg38 for bed file format
	#
	echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_coverage_graph_simplied.R $file $chrom_id 1 $size 'hg38' " | msub -q high -l nodes=1:ppn=1,mem=16gb -V -d $path -N graph$id -o out$id -e err$id -A proj-dm0001
	#echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_high_cov_graph.R $file $chrom_id 1 $size 'hg38' " | msub -q high -l nodes=1:ppn=1,mem=16gb -V -d $path -N graph$id -o out$id -e err$id -A proj-dm0001

	id=$((id+1))
done

#for file in `ls $path/BLI089*WGS*between*x_REPORT.txt`
for file in `ls $path/*WGS*between*_REPORT.txt`
#for file in `ls $path/*between1x_150x_REPORT.txt`
do
	# for hg37 for bed file format
	#
	#echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_coverage_graph_simplied.R $file $chrom_id 1 $size 'hg37' " | msub -q high -l nodes=1:ppn=1,mem=16gb -V -d $path -N graph$id -o out$id -e err$id -A proj-dm0001
	#echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_high_cov_graph.R $file $chrom_id 1 $size 'hg38' " | msub -q high -l nodes=1:ppn=1,mem=16gb -V -d $path -N graph$id -o out$id -e err$id -A proj-dm0001

	id=$((id+1))
done
