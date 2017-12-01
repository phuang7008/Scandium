#!/bin/bash

# This script is used to setup the eMerge annotation database for sequencing analysis
# it will use the annotation from the following directory:
# /hgsc_software/production/users/cbuhay/ExCID/working/ExCID_v2.1.5/database

# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/annotation_from_java_excid/setup_eMerge_annotation_databases.sh target_bed_file output_dir"
	exit
fi

# Capture target bed file and get rid of the 'chr' in front chromosome id
target_bed_file=$1

# get the working directory and cd to the directory
BASEDIR=$2
cd $BASEDIR
pwd
#echo `$BASEDIR`
printf "$BASEDIR\n"

# Before we continue, we need to sort the target_bed_file
#
target_base_name=`basename ${target_bed_file}`
echo "target base name is $target_base_name"
target_sorted_bed_file=$BASEDIR"/$target_base_name"_sorted
echo "sort $target_bed_file to produce $target_sorted_bed_file"
#sort_flag=" -k1,1 -V -s -k2,2n -k3,3n "
sort_flag=" -k1,1 -V -s -k4,4 -k2,2n "
sort $sort_flag -o $target_sorted_bed_file $target_bed_file

# For CDS database, dump everything into MySQL database named: eMerge_CDS37
#
echo "dump all exons/cds into MySQL database eMerge_CDS37 ==> eMergeCDSCoords.pl $target_sorted_bed_file"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/annotation_from_java_excid/eMergeCDS_SNPCoords.pl $target_sorted_bed_file "hg37"
#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/annotation_from_java_excid/eMergeCDSCoords.pl $target_sorted_bed_file "hg37"

# Now handle the real annotation database
# here, we need to sort the original file in a different order
# and then merged those CDSs with perfect coordinates
#
#sort_flag=" -k1,1 -V -s -k2,2n "
#target_merged_bed_file=$BASEDIR"/$target_base_name"_merged
#sort $sort_flag $target_bed_file | /stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 - | uniq - > $target_merged_bed_file

# now dump the merged results into the eMerge_Annotation37 database
#
#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/annotation_from_java_excid/eMergeAnnotation.pl $target_sorted_bed_file "hg37"

####
#END
####
