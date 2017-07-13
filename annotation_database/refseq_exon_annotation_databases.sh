#!/bin/bash

# This script is used to setup the gene intronic annotation database for sequencing analysis
# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/refseq_exon_annotation_databases.sh annotation_file target_bed_file"
	exit
fi

# the input file (BED_FILES) should contains the full path to the RefSeq, CCDS and 
# GenCode bed files downloaded from UCSC website
BED_FILES=$1

# Capture target bed file
target_bed_file=$2

# need to combine all the bed files into a single bed file
#
refseq_bed_file=""

while IFS='' read -r line || [[ -n "$line" ]];
do
	echo "$line"
	renamed_file=`basename $line`_renamed

	if [[ "${line,,}" == *"refseq"* ]]; then
		# here we are only interested in refseq source of annotation.
		# since some refseq names are the same but at different chromosome locations, we need to separate them.
		# For example: The name will be changed from NM_000015 to NM_000015-1, NM_000015-2 etc
		#
		/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/renameRefSeq.pl "$line" > "$renamed_file"

		refseq_bed_file="$renamed_file"_bed
		/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/generate_refseq_exon_bed.py -i "$renamed_file" > "$refseq_bed_file"
	else
		continue
	fi

done < "$BED_FILES"

# Before we continue, we need to sort the refseq_bed_files 
#
refseq_sorted_bed_file="$refseq_bed_file"_sorted
echo "sort $refseq_sorted_bed_file to produce $refseq_sorted_bed_file"
sort_flag="-k1,1 -k2,2n -k3,3n"
sort $sort_flag -o $refseq_sorted_bed_file $refseq_bed_file

# Next, we need to merge those exons with perfect coordinates.
#
refseq_sorted_bed_file_merged_perfect_matches="$refseq_sorted_bed_file"_merged_perfect_matches
echo "merge perfect matched regions in $refseq_sorted_bed_file to produce $refseq_sorted_bed_file_merged_perfect_matches"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$refseq_sorted_bed_file" | uniq - > "$refseq_sorted_bed_file_merged_perfect_matches"

# Here we are going to Get the intersect regions between exons and targets
# ====> this won't work at this moment as chromsome id is named differently between hg37 and hg38, need to find a way to handle this
exon_target_intersect_file="exon_target_intersect_for_gene_percentage_annotation"
echo "bedtools intersect to produce $exon_target_intersect_file from $refseq_sorted_bed_file_merged_perfect_matches and $target_bed_file"
/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $refseq_sorted_bed_file_merged_perfect_matches -b $target_bed_file > $exon_target_intersect_file

# Finally, dump everything into MySQL database named: Gene_RefSeq_Exon
#
echo "dump all exons into MySQL database Gene_RefSeq_Exon38 ==> geneExonAnnotation.pl $exon_target_intersect_file"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/geneExonAnnotation.pl "$exon_target_intersect_file"

####
#END
####




