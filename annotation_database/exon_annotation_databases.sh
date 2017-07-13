#!/bin/bash

# This script is used to setup the gene/exon annotation database for sequencing analysis
# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/exon_annotation_databases.sh annotation_file hgnc_version_hg38"
	exit
fi

# the input file (BED_FILES) should contains the full path to the RefSeq, CCDS and 
# GenCode bed files downloaded from UCSC website
BED_FILES=$1

# for HGNC official gene_symbol and dump them into the MySQL database
#
HGNC=$2
echo "processHGNCtoDB.pl $HGNC"
#/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/processHGNCtoDB.pl "$HGNC"

# need to combine all the bed files into a single bed file
#
combined_bed_files="all_exons_bed"

while IFS='' read -r line || [[ -n "$line" ]];
do
	echo "$line"
	tmp_file=`basename $line`_bed

	#if [[ $line == *"miRNA"* ]]; then
	if [[ "${line,,}" == *"mirna"* ]]; then
		/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/rearrange_miRNA_bed.py -i "$line" > "$tmp_file"
	else
		/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/generate_bed_file.py -i "$line" > "$tmp_file"
	fi

	# combine all the bed files together
	`cat $tmp_file >> $combined_bed_files`
done < "$BED_FILES"

# Before we continue, we need to sort the combined_bed_files 
#
combined_sorted_bed_file="$combined_bed_files"_sorted
echo "sort $combined_sorted_bed_file to produce $combined_sorted_bed_file"
sort_flag="-k1,1 -k2,2n -k3,3n"
sort $sort_flag -o $combined_sorted_bed_file $combined_bed_files

# Next, we need to merge those exons with perfect coordinates from different annotation sources.
#
combined_sorted_bed_file_merged_perfect_matches="$combined_sorted_bed_file"_merged_perfect_matches
echo "merge perfect matched regions in $combined_sorted_bed_file to produce $combined_sorted_bed_file_merged_perfect_matches"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$combined_sorted_bed_file" | uniq - > "$combined_sorted_bed_file_merged_perfect_matches"

# Now it is the time to generate exon partition file
#
exons_partitioned="$combined_sorted_bed_file"_merged_and_partitioned
echo "partition $combined_sorted_bed_file_merged_perfect_matches to produce $exons_partitioned using bedops"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedops -p "$combined_sorted_bed_file_merged_perfect_matches" | /stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap	--echo	--echo-map-id --delim '\t' - $combined_sorted_bed_file_merged_perfect_matches | uniq - > "$exons_partitioned"

# Finally, dump everything into MySQL database named: Exon_Regions38
#
echo "dump all partitioned exons into MySQL database Exon_Regions38 ==> exonAnnotations.pl $exons_partitioned"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/exonAnnotations.pl "$exons_partitioned"

####
#END
####




