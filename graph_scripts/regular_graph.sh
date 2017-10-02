#!/bin/bash

# this script is used to generate the OMIM gene transcript bed file from OMIM database
#
if [[ $# -ne 6 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/scripts/regular_graph.sh input_output_directory file_suffix_to_be_processed chrom_id lower_bound higher_bound type"
	exit
fi

# get the working directory and cd to the directory
BASEDIR=$1
cd $BASEDIR
printf "$BASEDIR\n"

# file suffix that we are interested in
SUFFIX=$2

# chromosome id we will draw graph on
CHROM_ID=$3

# graph boundaries
LOW=$4
HIGH=$5
TYPE=$6

id=0

for file in `ls $BASEDIR/*$SUFFIX`
do
	id=$((id+1))

	#echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/dev/make_graphs/R/low_high_coverage_log.R $file '1' 0 20 'hg19' " | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $path -N graph -o out_cap -e err_cap -A proj-dm0001

	echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_coverage_graph.R $file $CHROM_ID $LOW $HIGH $TYPE " | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $BASEDIR -N graph$id -o out$id -e err$id -A proj-dm0001
done
