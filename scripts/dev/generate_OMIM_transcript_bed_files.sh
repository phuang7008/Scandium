#!/bin/bash

# this script is used to generate the OMIM gene transcript bed file from OMIM database
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/scripts/generate_OMIM_transcript_bed_files.sh input_output_directory OMIM_database_file"
	exit
fi

# get the working directory and cd to the directory
BASEDIR=$1
cd $BASEDIR
printf "$BASEDIR\n"

OMIM_database_file=$2

# create an OMIM transcript bed file from OMIM exon bed file
OMIM_Transcript_bed="OMIM_transcripts.bed"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/create_new_OMIM_bed_file.pl $OMIM_database_file $OMIM_Transcript_bed

# here I need to calculate the GC content for each OMIM gene
# and get columns that we are interested in: chrom_id ($1), start ($2), end ($3), annotation ($4 not needed), GC content ($6)
#
OMIM_transcript_profile="OMIM_Transcript_Profile"
/hgsc_software/BEDTools/latest/bin/nucBed -fi /stornext/snfs5/next-gen/Illumina/bwa_references/h/hs37d5/hs37d5.fa -bed $OMIM_Transcript_bed | awk '{OFS="\t"; print $1,$2,$3,$6}' | tail -n +2  > $OMIM_transcript_profile

# Before we continue, we need to sort the OMIM_Transcript_Profile file
#
OMIM_Transcript_bed_sorted=$OMIM_Transcript_bed"_sorted"
sort_flag="-k1,1 -V -s -k2,2n -k3,3n"
sort $sort_flag -o $OMIM_Transcript_bed_sorted $OMIM_transcript_profile

# for Low_GC content 
#
OMIM_Low_GC="OMIM_Low_GC_Transcript.bed"
awk -F"\t" '$4<0.35' $OMIM_Transcript_bed_sorted > $OMIM_Low_GC

# for High GC content
#
OMIM_High_GC="OMIM_High_GC_Transcript.bed"
awk -F"\t" '$4>0.65' $OMIM_Transcript_bed_sorted > $OMIM_High_GC

# for Ave GC content
#
OMIM_Ave_GC="OMIM_Ave_GC_Transcript.bed"
awk -F"\t" '$4<0.5 && $4>0.45' $OMIM_Transcript_bed_sorted > $OMIM_Ave_GC

####
#END
####
