#!/bin/bash

# this script is used to generate the gene coverage graphs using gvcf file format
# the gene list is from OMIM database
#
if [[ $# -ne 7 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/scripts/generate_OMIM_genes_cov_graphs.sh OMIM_bed_file chrom_id input_output_directory suffix_low_cov_report_file lower_bound higher_bound type"
	exit
fi

# get the working directory and cd to the directory
BASEDIR=$3
cd $BASEDIR
printf "$BASEDIR\n"

# file suffix that we are interested in
SUFFIX=$4

# chromosome id we will draw graph on
CHROM_ID=$2

# graph boundaries
LOW=$5
HIGH=$6
TYPE=$7

OMIM_Transcript_bed=$1

# loop through the file list with the same SUFFIX intersect with OMIM transcript sorted bed file
# we will combine all the positions from the chrom id user interested in to a single file
#
combined_low_cov_position_file="combined_low_cov_position_for_chrom_"$CHROM_ID

# now sort the combined_low_cov_position_file
sorted_combined_file=$combined_low_cov_position_file"_sorted"
#sort_flag="-k1,1 -V -s -k2,2"
#sort $sort_flag -o $sorted_combined_file $combined_low_cov_position_file

# now we need to merge all regions that overlap each other 1%
merged_outfile=$sorted_combined_file"_merged"
#/hgsc_software/BEDTools/latest/bin/bedtools merge -i "$sorted_combined_file" | uniq - > $merged_outfile

# get the total length info from combined_low_cov_position_file
#
max_length=`/stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/get_total_length_to_draw.pl $merged_outfile $CHROM_ID`

echo "the max value is $max_length"

# go through the loop again to draw the graph
#
id=1
for file in `ls $BASEDIR/*$SUFFIX`
do
	#outfile=$file"_OMIM_only_size_recalculated"
	infile=$file"_OMIM_only"

	# draw graph for CHROM_ID
	echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/OMIM_gene_coverage_graph.R $merged_outfile $infile '1' $LOW $HIGH $max_length 1" | msub -q normal -l nodes=1:ppn=1,mem=15gb -V -d $BASEDIR -N graph -o out_cap$id -e err_cap$id -A proj-dm0001

	# echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/gene_coverage_graph.R $infile '1' $LOW $HIGH $max_length 1" | msub -q normal -l nodes=1:ppn=1,mem=12gb -V -d $BASEDIR -N graph -o out_cap$id -e err_cap$id -A proj-dm0001

	#echo "/hgsc_software/R/R-3.2.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/low_high_coverage_interval_HPC_by_name.R $outfile '1' $LOW $HIGH $LENGTH $TYPE" | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $BASEDIR -N graph -o out_cap$id -e err_cap$id -A proj-dm0001

	id=$((id+1))
done

####
#END
####
