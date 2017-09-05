#!/bin/bash

# This script is used to setup the gene/exon annotation database for sequencing analysis
# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/setup_exon_annotation_databases.sh miRNA_file output_directory"
	exit
fi

# get the working directory and cd to the directory
BASEDIR=$2
cd $BASEDIR
pwd
#echo `$BASEDIR`
printf "$BASEDIR\n"

# since miRNA file can't be downloaded using rsync or ftp, we have to download it first and tell the script where to find it
# things might change in the future and we might be able to use rsync or ftp in the future! 
# But for now, we need to do it manually!
miRNA_FILE=$1

# now process miRNA file
new_miRNA_file=`basename ${miRNA_FILE}`_simplified.bed
awk -F "\t| " '{print $1"\t"$4"\t"$5"\t"$9}' $miRNA_FILE | grep 'chr' | tr -s ';' '\t' | cut -f1,2,3,4,6 | sed s/ID=// | sed s/Name=// | sed s/chr//gi | awk '{t=$5; $5=$4; $4=t; print}' | awk '{ $4=$4"="$5; print}' > $new_miRNA_file
printf "Producing simplified miRNA file $new_miRNA_file \n"

# for HGNC official gene_symbol and dump them into the MySQL database
#
rsync -a -P rsync://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt .
HGNC='hgnc_complete_set.txt'
HGNC_simplified='hgnc_simplified'
cat $HGNC | cut -f 1,2,9,11,20,22,24,25,33 | sed s/HGNC:// | tr -d '"' | sed s/\|/,/g > $HGNC_simplified
printf "dump $HGNC_simplified info into DB using processHGNCtoDB.pl \n"
#/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/processHGNCtoDB.pl "$HGNC_simplified"

# now get rid of .txt file extension before preceed, as we will unzip the tar file into .txt
for f in $BASEDIR/*.txt; do mv -- "$f" "${f%.txt}"; done

#if [ ${HGNC: -4} == ".txt" ]; then
#	mv $HGNC "hgnc_complete_set"
#fi

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
    ls $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/generate_bed_file.py -i "$FILE.tmp" > "$FILE.bed"; done;
elif [ "$(uname)" == "Linux" ]; then
    ls --color=never $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls --color=never $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/generate_bed_file.py -i "$FILE.tmp" > "$FILE.bed"; done;
fi

# remove all the tmp files files
for f in $BASEDIR/*.tmp; do rm -f $f; done
#for f in $BASEDIR/*.txt; do rm -f $f; done

# save all the .txt.gz downloaded source files in case we need it in the future!
tmp_timestamp=$(date +%Y'_'%m'_'%d)
printf "time stamp is $tmp_timestamp\n"
for f in $BASEDIR/*.txt; do gzip $f; mv "$f.gz" "$f.$tmp_timestamp.gz"; done
mv "hgnc_complete_set" "hgnc_complete_set.txt"

# need to combine all the bed files into a single bed file and then remove the original bed files
#
combined_bed_files="all_exons_bed"
for f in $BASEDIR/*.bed; do cat $f >> $combined_bed_files; rm -f $f; done

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
printf "partition $combined_sorted_bed_file_merged_perfect_matches to produce $exons_partitioned using bedops \n"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedops -p "$combined_sorted_bed_file_merged_perfect_matches" | /stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap	--echo	--echo-map-id --delim '\t' - $combined_sorted_bed_file_merged_perfect_matches | uniq - > "$exons_partitioned"

# Finally, dump everything into MySQL database named: Exon_Regions38
#
printf "dump all partitioned exons into MySQL database Exon_Regions38 ==> exonAnnotations.pl $exons_partitioned \n"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/exonAnnotations.pl "$exons_partitioned"

####
#END
####
