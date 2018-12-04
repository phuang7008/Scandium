#!/bin/bash

# this is a wrapper script that is used to generate a whole chromosome coverage graph for viewing
# In order to run this script, you need to run the C/C++ ExCID first, with -U and -H options specified
# The best range should be between 1-100 (tested already). ie, -U 100 -H 1 (yes, NOT -L option)
#
# This will generate a file called XXX_WGS_between1x_100x_REPORT.txt or XXX_Capture_between1x_100x_REPORT.txt
# These generated report files will be in bed format and should be the one used for this script input files
# You don't have to specify them one by one, you could just specify the suffix of these files
#
if [[ $# -ne 7 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/regular_graph.sh input_directory output_directory file_suffix_to_be_processed(generated from C/C++ ExCID) chrom_id lower_bound(larger than 0)  higher_bound type"
	exit
fi

# get the working input directory and cd to the directory
INPUTDIR=$1

# this is the working output directory
BASEDIR=$2
cd $BASEDIR
printf "output directory is $BASEDIR\n"

# file suffix that we are interested in
SUFFIX=$3

# chromosome id we will draw graph on
CHROM_ID=$4

# graph boundaries
LOW=$5
HIGH=$6
TYPE=$7

id=0

for file in `ls $INPUTDIR/*$SUFFIX`
do
	id=$((id+1))

	#echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/dev/make_graphs/R/low_high_coverage_log.R $file '1' 0 20 'hg19' " | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $path -N graph -o out_cap -e err_cap -A proj-dm0001

	echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_coverage_graph.R $file $CHROM_ID $LOW $HIGH $TYPE " | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $BASEDIR -N graph$id -o out$id -e err$id -A proj-dm0001
done
