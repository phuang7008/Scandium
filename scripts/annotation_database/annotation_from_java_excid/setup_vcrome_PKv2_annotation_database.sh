#!/bin/bash

# This script is used to setup the VCRome+PKv2 annotation database for sequencing analysis
# it will use the annotation from the following directory:
# /hgsc_software/production/users/cbuhay/ExCID/working/ExCID_v2.1.5/database

# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 6 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/annotation_from_java_excid/setup_vcrome_PKv2_annotation_database.sh miRAN_annotation CCDS_annotation RefSeq_annotation VEGA_annotation PKv2_Coords output_dir"
	exit
fi

# Capture target bed file and get rid of the 'chr' in front chromosome id
miRNA_bed_file=$1
CCDS_bed_file=$2
Refseq_bed_file=$3
VEGA_bed_file=$4

# for PKv2 Coordinates
#
target_bed_file=$5

# get the working directory and cd to the directory
BASEDIR=$6
cd $BASEDIR
pwd
#echo `$BASEDIR`
printf "$BASEDIR\n"

# we are going to process them one at a time
# Before we continue, we need to sort them first
#
combined_bed_file="combined_bed"

for f in $miRNA_bed_file $CCDS_bed_file $Refseq_bed_file $VEGA_bed_file
do
	target_base_name=`basename ${f}`
	echo "target base name is $target_base_name"
	reformat_bed_file=$BASEDIR"/$target_base_name"_reformatted
	echo "Do awk to reformat the file $target_base_name"
	awk '{OFS="\t"; print$1,$2,$3,$4"="$5"="$6}' $f > $reformat_bed_file

	# need to combine all the bed files
	#
	cat $reformat_bed_file >> $combined_bed_file
done

# Here we need to re-sort everything in the combined_bed_file
#
sort_flag=" -k1,1 -V -s -k2,2n -k3,3n "
combined_bed_sorted="combined_bed_sorted"
echo "Sort $combined_bed_file to $combined_bed_sorted"
sort $sort_flag -o $combined_bed_sorted $combined_bed_file

# next merged the perfect matched ones
#
combined_bed_merged="combined_bed_sorted_merged"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$combined_bed_sorted" | uniq - > $combined_bed_merged

# Now it is the time to generate exon partition file
#
exons_partitioned="combined_bed_sorted_merged_and_partitioned"
printf "partition $combined_bed_merged to produce $exons_partitioned using bedops (only works for chr1-9 and chrX-Y\n"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedops -p "$combined_bed_merged" | /stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo  --echo-map-id --delim '\t' - $combined_bed_merged | uniq - > "$exons_partitioned"

# BIG NOTE:
# It seems that the above command doesn't work for chromosome 10-22, thus, I need to do extra things here to handle this
#
echo "extract everything from 1-9 and X and Y that have been annotated"
worked_chromosome_partitioned_file="annotation_for_1_9_X_Y"
awk -F"\t" '$1<10 || $1>22' $exons_partitioned> $worked_chromosome_partitioned_file

echo "extract everything from 10-22 that haven't been annotated"
partitioned_chr10_22_only="chr10_22_partitioned_bed"
awk -F"\t" '$1>9' $exons_partitioned  | awk -F"\t" '$1<=22' | cut -f1,2,3 > $partitioned_chr10_22_only

echo "get annotations from merged file"
merged_annotation_chr10_22_only="merged_annotation_chr10_22_only"
awk -F"\t" '$1>9' $combined_bed_merged | awk -F"\t" '$1<=22' > $merged_annotation_chr10_22_only

echo "now do the bedmap to generated annotated info for chr10 to chr22"
annotation_for_chr10_22="annotation_for_chr10_22"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo --echo-map-id --delim '\t' $partitioned_chr10_22_only $merged_annotation_chr10_22_only > $annotation_for_chr10_22

echo "combine them together!"
final_partitioned_file="final_partitioned_annotations.bed"
cat $worked_chromosome_partitioned_file > $final_partitioned_file
cat $annotation_for_chr10_22 >> $final_partitioned_file

# Now, dump everything into MySQL database named: VCRomePKv2_Annotation37
#
echo "dump all exons/cds into MySQL database VCRomePKv2_Annotation37/38	==> vcromePKv2Annotation.pl $target_sorted_bed_file 'hg37' $counter"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/annotation_from_java_excid/vcromePKv2Annotation.pl $final_partitioned_file "hg37"

##################################################################

# Now we need to process the genes' CDS coordinates
# Note: for VCRome + PKv2, we are only interested in the coordinates for PKv2 part, not VCRome part!
# Before we continue, we need to sort the target_bed_file
#
target_base_name=`basename ${target_bed_file}`
echo "target base name is $target_base_name"
target_sorted_bed_file=$BASEDIR"/$target_base_name"_sorted
echo "sort $target_bed_file to produce $target_sorted_bed_file"
sort_flag=" -k1,1 -V -s -k4,4 -k2,2n "
sort $sort_flag -o $target_sorted_bed_file $target_bed_file
echo "dump all exons/cds into MySQL database PKv2_CDS37 ==> PKv2CDSCoords.pl $target_sorted_bed_file"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/annotation_from_java_excid/PKv2CDSCoords.pl $target_sorted_bed_file "hg37"

####
#END
####
