#!/bin/bash

# This script is used to setup the eMerge annotation database for sequencing analysis
# it will use the annotation from the following directory:
# /hgsc_software/production/users/cbuhay/ExCID/working/ExCID_v2.1.5/database

# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 3 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/annotation_from_java_excid/setup_right10k_annotation_databases.sh gene_bed_file snp_bed_file output_dir"
	exit
fi

# Capture target bed file and get rid of the 'chr' in front chromosome id
gene_bed_file=$1
snp_bed_file=$2

# get the working directory and cd to the directory
BASEDIR=$3
cd $BASEDIR
pwd
#echo `$BASEDIR`
printf "$BASEDIR\n"

# we are going to process them one at a time
# Before we continue, we need to sort them first
#
target_base_name=`basename ${gene_bed_file}`
echo "target base name is $target_base_name"
target_sorted_bed_file=$BASEDIR"/$target_base_name"_sorted
echo "sort $target_bed_file to produce $target_sorted_bed_file"
#sort_flag=" -k1,1 -V -s -k2,2n -k3,3n "
sort_flag=" -k1,1 -V -s -k4,4 -k2,2n "
sort $sort_flag -o $target_sorted_bed_file $gene_bed_file

# Now, dump everything into MySQL database named: Right10K_CDS37
#
echo "dump all exons/cds into MySQL database Right10K_CDS37 ==> right10kCDSCoords.pl $target_sorted_bed_file 'hg37' 'gene'"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/annotation_from_java_excid/right10kCDSCoords.pl $target_sorted_bed_file "hg37" "gene"

# for SNP
snp_base_name=`basename ${snp_bed_file}`
echo "snp basename is $snp_base_name"
snp_sorted_bed_file=$BASEDIR"/$snp_base_name"_sorted
echo "sort $snp_base_name to produce $snp_sorted_bed_file"
sort $sort_flag -o $snp_sorted_bed_file $snp_bed_file

# dump snp into database
echo "dump all SNPs into MySQL database Right10K_SNP37 ===> right10kCDSCoords.pl $snp_sorted_bed_file 'hg37' 'snp'" 
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/annotation_from_java_excid/right10kCDSCoords.pl $snp_sorted_bed_file "hg37" "snp"

####
#END
####
