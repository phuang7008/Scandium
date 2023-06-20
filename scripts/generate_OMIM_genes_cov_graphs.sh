#!/bin/bash

# this script is used to generate the gene coverage graphs. 
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

# create an OMIM transcript bed file from OMIM exon bed file
#OMIM_Transcript_bed="OMIM_transcripts.bed"
#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/create_new_OMIM_bed_file.pl $OMIM_Transcript_bed

# Before we continue, we need to sort the OMIM_Transcript_bed 
# this will be pre-processed!
#
OMIM_Transcript_bed_sorted=$OMIM_Transcript_bed"_sorted"
sort_flag="-k1,1 -V -s -k2,2"
sort $sort_flag -o $OMIM_Transcript_bed_sorted $OMIM_Transcript_bed

# loop through the file list with the same SUFFIX intersect with OMIM transcript sorted bed file
# and get the max length from all files
#
max_length=1000

for file in `ls $BASEDIR/*$SUFFIX`
do
	# need to remove first 3 lines in the SUFFIX files, as they are comments only
	# Note: use this only when you are using the results from new EXCID coverage report
	# for other files, if there is no comment lines, don't use this
	#
	Comments_Removed=$file"_tmp"
	#tail -n +4 $file > $Comments_Removed
	sed '$d' $file > $Comments_Removed

	# now need to do the intersect
	outfile=$file"_OMIM_only"
	/hgsc_software/BEDTools/latest/bin/bedtools intersect -b $Comments_Removed -a $OMIM_Transcript_bed_sorted -wb > $outfile
	#/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $Comments_Removed -b $OMIM_Transcript_bed_sorted > $outfile

	# after intersect, we need to add coverage information onto the annotation field (field #4)
	#
	outfile_cov=$file"_OMIN_w_Cov"
	`awk '{OFS="\t"; print $1,$2,$3,$4"-"$8}' $outfile > $outfile_cov`

	# after intersect, we need to merge those with exact same coordinates from different transcripts
	#
	#merged_outfile=$outfile"_merged"
	#/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$outfile" | uniq - > $merged_outfile

	# after intersect/merge, the sizes of low/high coverage regions will be changed, here we will re-calculate them accordingly
	#
	#recalculated_file=$merged_outfile"_size_recalculated"
	recalculated_file=$outfile_cov"_size_recalculated"
	`sed  -e "s/\-/\\t/g" $outfile_cov | awk '{OFS="\t"; print $1,$2,$3,$3-$2+1,$8,$4,$5,$6,$7}'  > $recalculated_file`

	# get the total length info
	LENGTH=`/stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/get_total_length_to_draw.pl $recalculated_file $CHROM_ID`
	echo "The length is $LENGTH"

	if [ "$max_length" -lt "$LENGTH" ]; then
		#echo "alignment is on! $LENGTH"
		max_length=$LENGTH
	fi
done

echo "the max value is $max_length"

# go through the loop again to draw the graph
#
id=1
for file in `ls $BASEDIR/*$SUFFIX`
do
	#outfile=$file"_OMIM_only_size_recalculated"
	infile=$file"_OMIN_w_Cov_size_recalculated"

	# draw graph for CHROM_ID
	echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/gene_coverage_graph.R $infile '1' $LOW $HIGH $max_length 1" | msub -q normal -l nodes=1:ppn=1,mem=12gb -V -d $BASEDIR -N graph -o out_cap$id -e err_cap$id -A proj-dm0001

	#echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/low_high_coverage_interval_HPC_by_name.R $outfile '1' $LOW $HIGH $LENGTH $TYPE" | msub -q normal -l nodes=1:ppn=1,mem=32gb -V -d $BASEDIR -N graph -o out_cap$id -e err_cap$id -A proj-dm0001

	id=$((id+1))
done

####
#END
####
