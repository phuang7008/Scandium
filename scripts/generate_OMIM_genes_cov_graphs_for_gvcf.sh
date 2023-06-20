#!/bin/bash

# this script is used to generate the gene coverage graphs using gvcf file format
# the gene list is from OMIM database
#
if [[ $# -ne 8 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/scripts/generate_OMIM_genes_cov_graphs.sh OMIM_bed_file chrom_id input_directory output_directory suffix_low_cov_report_file lower_bound higher_bound DB_version"
	exit
fi

# get the working input and output directories and cd to the directory
INPUTDIR=$3
BASEDIR=$4		# output directory
cd $BASEDIR
printf "$BASEDIR\n"

# file suffix that we are interested in
SUFFIX=$5

# chromosome id we will draw graph on
CHROM_ID=$2

# graph boundaries
LOW=$6
HIGH=$7
DB_Version=$8

OMIM_Transcript_bed=$1

# set the sort flag
#
sort_flag="-k1,1 -V -s -k2,2n -k3,3n"

# loop through the file list with the same SUFFIX and intersect with OMIM transcript sorted bed file
# we will combine all the positions from the chrom id user interested in to a single file
#
combined_low_cov_position_file="combined_low_cov_position_for_chrom_"$CHROM_ID
longest=0

for file in `ls $INPUTDIR/*$SUFFIX`
do
	# now need to do the intersect
	#
	outfile=$file"_OMIM_only"
	/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $file -b $OMIM_Transcript_bed | awk -v a="$CHROM_ID" '{OFS="\t"; if ($1==a) print}' | sort $sort_flag | uniq > $outfile

	tmp_longest=`less $outfile | wc -l `
	if [ $tmp_longest -gt $longest ]; then
		longest=$tmp_longest
	fi

	# we can't merge it as they have different coverage number, unless we use average for combined rows
	# now merge it, with gap of 1 (such as 500 and 501 will be merged)
	# 
	merged_file=$file"_OMIM_merged"
	/hgsc_software/BEDTools/latest/bin/bedtools merge -i $outfile -d 1 -c 4 -o mean > $merged_file

	# dump current library low coverage regions with current chromosome id to a combined low coverage region file
	# the combined low coverage regions will form the backbone of our graph
    #
    awk -v a="$CHROM_ID" '{OFS="\t"; if ($1==a) print $1,$2,$3}' $outfile >> $combined_low_cov_position_file

	echo "The total low coverage regions for current library $file is $longest === DONE!"
done

# now sort the combined_low_cov_position_file
#
sorted_combined_file=$combined_low_cov_position_file"_sorted"
sort $sort_flag -o $sorted_combined_file $combined_low_cov_position_file

# now we need to merge all regions that overlap each other (with gap of 1, such as 500 and 502 will be merged)
#
merged_combined_file=$sorted_combined_file"_merged"
/hgsc_software/BEDTools/latest/bin/bedtools merge -i "$sorted_combined_file" -d 2 | uniq - > $merged_combined_file

# get the total length info from combined_low_cov_position_file
#
max_length=`/stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/get_total_length_to_draw.pl $merged_combined_file $CHROM_ID`

echo "the max value is $max_length"

# go through the loop again to draw the graph
#
id=11
for file in `ls $INPUTDIR/*$SUFFIX`
do
	#infile=$file"_OMIM_only"
	merged_file=$file"_OMIM_merged"

	# modified Zhuoyi's suggestion
    # now create a new bed file using Perl and insert 0 as the coverage into all other positions on the entire chromosome
    #
    cov_on_whole_chrom=$file"_OMIM_gene_on_whole_chromosome_"$CHROM_ID
    /hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/pad_all_zeros_for_whole_chromosome.pl $merged_file $cov_on_whole_chrom $CHROM_ID $DB_Version

	# intersect with combined sorted merged file
	# 
	final_input_file=$file"_OMIM_alinged_across_all"
	/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $cov_on_whole_chrom -b $merged_combined_file > $final_input_file

	# draw graph for CHROM_ID, scaled to the union size of all (low or high) coverage positions
	#
	echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/OMIM_gene_coverage_union_graph_aligned_all.R  $final_input_file $CHROM_ID $LOW $HIGH $max_length 1" | msub -q normal -l nodes=1:ppn=1,mem=15gb -V -d $BASEDIR -N graph$id -o out_cap$id -e err_cap$id -A proj-dm0001

	# now draw draw for CHROM_ID, onto the entire chromosome
	#
	echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_coverage_graph.R $final_input_file $CHROM_ID $LOW $HIGH $DB_Version " | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $BASEDIR -N pics_W_$id -o out_W_$id -e err_W_$id -A proj-dm0001

	if [ $HIGH -ne 100 ]; then
		#echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_coverage_graph.R $infile $CHROM_ID $LOW 100 $DB_Version " | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $BASEDIR -N graph_100_$id -o out_100_$id -e err_100_$id -A proj-dm0001
		#echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/OMIM_gene_coverage_union_graph.R $merged_combined_file $infile $CHROM_ID $LOW 100 $max_length 1" | msub -q normal -l nodes=1:ppn=1,mem=15gb -V -d $BASEDIR -N graph_100_$id -o out_100_$id -e err_100_$id -A proj-dm0001
		echo "100x"
	fi

	if [ $HIGH -ne 25 ]; then
		#echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/regular_whole_chrom_coverage_graph.R $infile $CHROM_ID $LOW 25 $DB_Version " | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $BASEDIR -N pics_25_$id -o out_25_$id -e err_25_$id -A proj-dm0001
		#echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/OMIM_gene_coverage_union_graph.R $merged_combined_file $infile $CHROM_ID $LOW 25 $max_length 1" | msub -q normal -l nodes=1:ppn=1,mem=15gb -V -d $BASEDIR -N graph_25_$id -o out_25_$id -e err_25_$id -A proj-dm0001
		echo "25x"
	fi

	id=$((id+1))
done

####
#END
####
