#!/bin/bash

# this is a wrapper script that is used to generate a whole chromosome or Capture coverage graph for viewing
# In order to run this script, you need to run the C/C++ ExCID first. For Capture, you need to specify -t
# The ExCID run will generate XXX_WGS_cov.fasta or XXX_Capture_cov.fasta files
# The graphs generated will be at the base level resolution!
#
if [[ $# -ne 6 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/create_graph_from_cov_fasta_files.sh input_directory output_directory chrom_id lower_bound(larger than 0) higher_bound genome_version"
	exit
fi

# get the working input directory and cd to the directory
INPUTDIR=$1

# this is the working output directory
BASEDIR=$2
cd $BASEDIR
printf "output directory is $BASEDIR\n"

# chromosome id we will draw graph on
CHROM_ID=$3

# graph boundaries
LOW=$4
HIGH=$5
VERSION=$6

# for WGS
#
id=0

for file in `ls $INPUTDIR/*WGS_cov.fasta`
do
	id=$((id+1))

	echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/coverage_plot_from_cov_fasta_files.R $file $CHROM_ID $LOW $HIGH $VERSION 'wgs'" | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $BASEDIR -N graph$id -o out$id -e err$id -A proj-dm0001
done

# for Capture
#
for file in `ls $INPUTDIR/*Capture_cov.fasta`
do
	id=$((id+1))

	echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/coverage_plot_from_cov_fasta_files.R $file $CHROM_ID $LOW $HIGH $VERSION 'capture'" | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $BASEDIR -N graph$id -o out$id -e err$id -A proj-dm0001
done
