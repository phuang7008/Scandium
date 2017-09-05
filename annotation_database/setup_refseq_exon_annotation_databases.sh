#!/bin/bash

# This script is used to setup the gene-refseq-exon annotation database for sequencing analysis
# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/setup_refseq_exon_annotation_databases.sh target_bed_file output_dir"
	exit
fi

# Capture target bed file and get rid of the 'chr' in front chromosome id
original_target_bed_file=$1
target_bed_file=`basename $original_target_bed_file`.modified
cat $original_target_bed_file | sed 's/chr//g' > $target_bed_file

# get the working directory and cd to the directory
BASEDIR=$2
cd $BASEDIR
pwd
#echo `$BASEDIR`
printf "$BASEDIR\n"

# here we are only interested in refseq source of annotation.
# since some refseq names are the same but at different chromosome locations, we need to separate them.
# For example: The name will be changed from NM_000015 to NM_000015-1, NM_000015-2 etc
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz .

# untar the zipped files
renamed_file='RefSeq_renamed'
if [ "$(uname)" == "Darwin" ]; then
    ls $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/renameRefSeq.pl "$FILE.tmp" > $renamed_file; done;
elif [ "$(uname)" == "Linux" ]; then
    ls --color=never $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls --color=never $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/renameRefSeq.pl "$FILE.tmp" > $renamed_file; done;
fi

refseq_bed_file="$renamed_file"_bed
/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/generate_refseq_exon_bed.py -i "$renamed_file" > "$refseq_bed_file"

# save all the .txt.gz downloaded source files in case we need it in the future!
tmp_timestamp=$(date +%Y'_'%m'_'%d)
printf "time stamp is $tmp_timestamp\n"
for f in $BASEDIR/*.txt; do gzip $f; mv "$f.gz" "$f.$tmp_timestamp.gz"; done

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
# ====> this won't work at this moment as chromsome id is named differently between hg37 and hg38, need to find a way to handle this (solved)
# in addition, we need liftover version of VCRome+PKv2 for hg38 (done: /users/ws144320/scratch/liftover/HG38_lom_vcrome2.1_with_PKv2.bed)
exon_target_intersect_file="exon_target_intersect_for_gene_percentage_annotation"
echo "bedtools intersect to produce $exon_target_intersect_file from $refseq_sorted_bed_file_merged_perfect_matches and $target_bed_file"
/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $refseq_sorted_bed_file_merged_perfect_matches -b $target_bed_file > $exon_target_intersect_file

# Finally, dump everything into MySQL database named: Gene_RefSeq_Exon
#
echo "dump all exons into MySQL database Gene_RefSeq_Exon38 ==> geneExonAnnotation_liftover.pl $exon_target_intersect_file"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/geneExonAnnotation.pl "$exon_target_intersect_file"

####
#END
####
