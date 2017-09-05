#!/bin/bash

# This script is used to setup the gene intronic annotation database for sequencing analysis
# First, we need to run the exon_annotation_databases.sh to create sorted_merged_partitioned_exon_file
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/setup_intron_annotation_databases.sh Sorted_Merged_Partitioned_Exon_FILE output_directory"
	exit
fi

# get the working directory and cd to the directory
BASEDIR=$2
cd $BASEDIR
pwd
#echo `$BASEDIR`
printf "$BASEDIR\n"

# since miRNAs don't contain introns, so we don't have to worry about them here
partitioned_Sorted_Merged_Exon_FILE=$1

# for HGNC official gene_symbol and dump them into the MySQL database
# This part should be already done by setup_exon_annotation_databases.sh script
#
#rsync -a -P rsync://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt .
#HGNC='hgnc_complete_set.txt'
#HGNC_simplified='hgnc_simplified'
#cat $HGNC | cut -f 1,2,9,11,20,22,24,25,33 | sed s/HGNC:// | tr -d '"' | sed s/\|/,/g > $HGNC_simplified
#printf "dump $HGNC_simplified info into DB using processHGNCtoDB.pl \n"
#/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/processHGNCtoDB.pl "$HGNC_simplified"

# now get rid of .txt file extension before preceed, as we will unzip the tar file into .txt
#for f in $BASEDIR/*.txt; do mv -- "$f" "${f%.txt}"; done
#mv "hgnc_complete_set" "hgnc_complete_set.txt"

# now we need to download the newest gene annotation from various sources (such as refseq, ccds and gencode)
#
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz .
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ccdsGene.txt.gz .
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV26.txt.gz .
#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeCompV26.txt.gz .
#rsync -a -P rsync://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3 .

# untar the zipped files
if [ "$(uname)" == "Darwin" ]; then
    ls $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/extractIntrons.py -i "$FILE.tmp" > "$FILE.bed"; done;
elif [ "$(uname)" == "Linux" ]; then
    ls --color=never $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls --color=never $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/extractIntrons.py -i "$FILE.tmp" > "$FILE.bed"; done;
fi

# remove all the tmp files files
#for f in $BASEDIR/*.tmp; do rm -f $f; done
#for f in $BASEDIR/*.txt; do rm -f $f; done

# save all the .txt.gz downloaded source files in case we need it in the future!
tmp_timestamp=$(date +%Y'_'%m'_'%d)
printf "time stamp is $tmp_timestamp\n"
for f in $BASEDIR/*.txt; do gzip $f; mv "$f.gz" "$f.$tmp_timestamp.gz"; done

# need to combine all the bed files into a single bed file and then remove the original bed files
#
combined_bed_files="all_introns_bed"
for f in $BASEDIR/*.bed; do cat $f >> $combined_bed_files; rm -f $f; done

# Before we continue, we need to sort the combined_bed_files 
#
combined_sorted_bed_file="$combined_bed_files"_sorted
echo "sort $combined_sorted_bed_file to produce $combined_sorted_bed_file"
sort_flag="-k1,1 -k2,2n -k3,3n"
sort $sort_flag -o $combined_sorted_bed_file $combined_bed_files

# Next, we need to merge those introns with perfect coordinates from different annotation sources.
#
combined_sorted_bed_file_merged_perfect_matches="$combined_sorted_bed_file"_merged_perfect_matches
echo "merge perfect matched regions in $combined_sorted_bed_file to produce $combined_sorted_bed_file_merged_perfect_matches"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$combined_sorted_bed_file" | uniq - > "$combined_sorted_bed_file_merged_perfect_matches"

# Now it is the time to generate intron partition file
#
introns_partitioned="$combined_sorted_bed_file"_merged_and_partitioned
echo "partition $combined_sorted_bed_file_merged_perfect_matches to produce $introns_partitioned using bedops"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedops -p "$combined_sorted_bed_file_merged_perfect_matches" | /stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap	--echo	--echo-map-id --delim '\t' - $combined_sorted_bed_file_merged_perfect_matches | uniq - > "$introns_partitioned"

# After that, we need to remove parts of exons that overlaps with intron regions
#
final_intron_regions="$combined_sorted_bed_file"_FINAL
echo "General final intronic regions with annotation from $introns_partitioned by extracting exons from $partitioned_Sorted_Merged_Exon_FILE to produce the $final_intron_regions"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedops -d $introns_partitioned $partitioned_Sorted_Merged_Exon_FILE | /stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo --echo-map-id --delim '\t' - $introns_partitioned | uniq - > $final_intron_regions

# Finally, dump everything into MySQL database named: Intron_Regions38
#
echo "dump all partitioned introns into MySQL database Intron_Regions38 ==> forIntronicRegions.pl $introns_partitioned"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/forIntronicRegions.pl "$final_intron_regions"

####
#END
####
