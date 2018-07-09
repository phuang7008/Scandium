#!/bin/bash

path=`pwd`
id=1
size=200
type=hg38

file=$path/$1

echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_coverage_graph_simplied.R $file '7' 1 $size $type " | msub -q high -l nodes=1:ppn=1,mem=32gb -V -d $path -N graph$id -o out$id -e err$id -A proj-dm0001

#echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_coverage_graph.R $file '7' 1 $size $type " | msub -q high -l nodes=1:ppn=1,mem=32gb -V -d $path -N graph$id -o out$id -e err$id -A proj-dm0001
