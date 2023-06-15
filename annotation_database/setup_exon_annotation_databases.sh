#!/bin/bash

# This script is used to setup the gene/exon annotation database for sequencing analysis
# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 5 ]]; then
	echo "Illegal Number of Parameters"
	echo "/SCRIPT_PATH/setup_exon_annotation_databases.sh output_directory Perl_Path Bedops_Path db_version(hg38 or hg37) annotation_source_type"
    echo "annotation_source_type: refseq or gencode or all (refseq+gencode+ccds+microRNA and vega(only in hg37)"
	exit
fi

# we don't want to hard code the path to the script. 
# we will use the following to get the bash script path
# this needs to be done before I change the directory
#
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# get the working directory and cd to the directory
BASEDIR=$1
cd $BASEDIR
pwd
printf "$BASEDIR\n"

# get perl interpret
PERL_PATH=$2

# get bedops path
BEDOPS_PATH=$3

# get the gene annotation version
gene_db_version=$4

# get annotation source
annotation_source_type=$5

# now process miRNA file
new_miRNA_file=hsa_miRNA_simplified.bed

# For microRNA
if [ $annotation_source_type == "all" ]; then
    if [ "$gene_db_version" == "hg38" ]; then
	    wget https://mirbase.org/ftp/CURRENT/genomes/hsa.gff3
	    #printf "Producing simplified miRNA file $new_miRNA_file \n"
    else
        wget https://mirbase.org/ftp/20/genomes/hsa.gff3    # version 20 is the last hg37 available, 21+ are all hg38
    fi

	awk -F "\t| " '{print $1"\t"$4"\t"$5"\t"$9}' hsa.gff3 | grep '^chr' | tr -s ';' '\t' | cut -f1,2,3,4,6 | sed s/ID=// | sed s/Name=// | sed s/^chr//gi | awk '{t=$5"\t"; $5=$4"\t"; $4=t"\t"; print}' | awk '{ $4=$4"|exon_1|miRNA_1="$5; print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > $new_miRNA_file
fi

# for HGNC official gene_symbol and dump them into the MySQL database
#
rsync -a -P rsync://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt .
HGNC='hgnc_complete_set.txt'
HGNC_simplified='hgnc_simplified'
cat $HGNC | cut -f 1,2,9,11,20,21,22,24,25,33 | sed s/HGNC:// | tr -d '"' | sed s/\|/,/g > $HGNC_simplified

# when you first run it, please un-comment out the line: processHGNCtoDB.pl
# As I am testing, I don't want to do it over and over
#
printf "dump $HGNC_simplified info into DB using processHGNCtoDB.pl \n"
$SCRIPT_PATH/processHGNCtoDB.pl "$HGNC_simplified" "$gene_db_version"

# now get rid of .txt file extension before preceed, as we will unzip the tar file into .txt
#
for f in $BASEDIR/*.txt; do mv -- "$f" "${f%.txt}"; done

# now we need to download the newest gene annotation from various sources (such as refseq, ccds and gencode)
#
if [ "$gene_db_version" == "hg38" ]; then
	echo "For hg38"
	if [ $annotation_source_type == "all" ] || [ $annotation_source_type == "refseq" ]; then
        rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz .
    fi

    if [ $annotation_source_type == "all" ] || [ $annotation_source_type == "gencode" ]; then
	    rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV26.txt.gz .
	    #rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeCompV26.txt.gz .
    fi

    if [ $annotation_source_type == "all" ]; then
	    rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ccdsGene.txt.gz .
    fi
else
	echo "For hg19"
	if [ $annotation_source_type == "all" ] || [ $annotation_source_type == "refseq" ]; then
	    rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz .
	    #rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeCompV19.txt.gz .
    fi
	if [ $annotation_source_type == "all" ] || [ $annotation_source_type == "gencode" ]; then
	    rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz .
    fi

    if [ $annotation_source_type == "all" ]; then
	    rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ccdsGene.txt.gz .
	    rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vegaGene.txt.gz .
    fi
fi

# untar the zipped files
# it seems we need to use the CDS start/end for the final calculation, so we will need to add them here
#
if [ "$(uname)" == "Darwin" ]; then
	echo "For Darwin"
    ls $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
	ls $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; $PERL_PATH $SCRIPT_PATH/generate_bed_file.pl "$FILE.tmp" "$gene_db_version" | sed s/^chr//gi | awk '$1!~/_/ {print}' > "$FILE.bed"; done;
elif [ "$(uname)" == "Linux" ]; then
	echo "For Linux"
    ls --color=never $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
	ls --color=never $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; $PERL_PATH $SCRIPT_PATH/generate_bed_file.pl "$FILE.tmp" "$gene_db_version" | sed s/^chr//gi | awk '$1!~/_/ {print}' > "$FILE.bed"; done;
fi

# remove all the tmp files files
#
#for f in $BASEDIR/*.tmp; do rm -f $f; done
#for f in $BASEDIR/*.txt; do rm -f $f; done

# save all the .txt.gz downloaded source files in case we need it in the future!
# especially when there is a new development and we want to compare the results using the same version of gene annotation
#
tmp_timestamp=$(date +%Y'_'%m'_'%d)
printf "time stamp is $tmp_timestamp\n"
for f in $BASEDIR/*.txt; do gzip $f; mv "$f.gz" "$f.$tmp_timestamp.gz"; done
#mv "hgnc_complete_set" "hgnc_complete_set.txt"

# need to combine all the bed files into a single bed file and then remove the original bed files
#
combined_bed_files="all_exons_bed"
for f in $BASEDIR/*.bed; do cat $f >> $combined_bed_files; done 

# here we set the sort flags as it will be used multiple times
#
sort_flag="-k1,1 -V -s -k2,2n -k3,3n"

# Now we will use the individual bed file to create a database on its own!
# These individual databases are needed for the Batch Analysis as it allows users to 
# choose which database annotations he/she wants
#
for f in $BASEDIR/*.bed 
do
	echo "Now processing $f"

	# for each bed file, we need to sort them first
	#
	echo "Sorting bed file"
	tmp_sorted_bed=$f."sorted"
	sort $sort_flag -o $tmp_sorted_bed $f

	# now merge all exons with the exact same coordinates from different transcripts
	#
	echo "Merged regions for sorted bed file"
	tmp_sorted_merged_bed=$tmp_sorted_bed."merged"
	$BEDOPS_PATH/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$tmp_sorted_bed" | uniq - > $tmp_sorted_merged_bed
done

# Before we continue, we need to sort the combined_bed_files 
#
combined_sorted_bed_file="$combined_bed_files"_sorted
echo "sort $combined_sorted_bed_file to produce $combined_sorted_bed_file"
sort $sort_flag -o $combined_sorted_bed_file $combined_bed_files

# Next, we need to merge those exons with perfect coordinates from different annotation sources.
#
combined_sorted_bed_file_merged_perfect_matches="$combined_sorted_bed_file"_merged_perfect_matches
echo "merge perfect matched regions in $combined_sorted_bed_file to produce $combined_sorted_bed_file_merged_perfect_matches"
echo "$BEDOPS_PATH/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 $combined_sorted_bed_file | uniq - > $combined_sorted_bed_file_merged_perfect_matches"
$BEDOPS_PATH/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$combined_sorted_bed_file" | uniq - > "$combined_sorted_bed_file_merged_perfect_matches"

# Now it is the time to generate exon partition file
#
exons_partitioned="$combined_sorted_bed_file"_merged_and_partitioned
printf "partition $combined_sorted_bed_file_merged_perfect_matches to produce $exons_partitioned using bedops \n"
echo "$BEDOPS_PATH/bedops -p $combined_sorted_bed_file_merged_perfect_matches | $BEDOPS_PATH/bedmap	--echo	--echo-map-id --delim '\t' - $combined_sorted_bed_file_merged_perfect_matches | uniq - > $exons_partitioned"
$BEDOPS_PATH/bedops -p "$combined_sorted_bed_file_merged_perfect_matches" | $BEDOPS_PATH/bedmap	--echo --echo-map-id --delim '\t' - $combined_sorted_bed_file_merged_perfect_matches | uniq - > "$exons_partitioned"

# BIG NOTE:
# It seems that the above command doesn't work for chromosome 10-22, thus, I need to do extra things here to handle this
#
echo "extract everything from 1-9 and X and Y that have been annotated"
worked_chromosome_partitioned_file="annotation_for_1_9_X_Y"
#echo "awk -F\"\t\" '$1<10 || $1>22' $exons_partitioned> $worked_chromosome_partitioned_file"
awk -F"\t" '$1<10 || $1>22' $exons_partitioned> $worked_chromosome_partitioned_file

echo "extract everything from 10-22 that haven't been annotated"
partitioned_chr10_22_only="chr10_22_partitioned_bed"
#echo "awk -F\"\t\" '$1>9' $exons_partitioned  | awk -F\"\t\" '$1<=22' | cut -f1,2,3 > $partitioned_chr10_22_only"
awk -F"\t" '$1>9' $exons_partitioned  | awk -F"\t" '$1<=22' | cut -f1,2,3 > $partitioned_chr10_22_only

echo "get annotations from merged file"
merged_annotation_chr10_22_only="merged_annotation_chr10_22_only"
#echo "awk -F\"\t\" '$1>9' $combined_sorted_bed_file_merged_perfect_matches | awk -F\"\t\" '$1<=22' > $merged_annotation_chr10_22_only"
awk -F"\t" '$1>9' $combined_sorted_bed_file_merged_perfect_matches | awk -F"\t" '$1<=22' > $merged_annotation_chr10_22_only

echo "now do the bedmap to generated annotated info for chr10 to chr22"
annotation_for_chr10_22="annotation_for_chr10_22"
#echo "$BEDOPS_PATH/bedmap --echo --echo-map-id --delim '\t' $partitioned_chr10_22_only $merged_annotation_chr10_22_only > $annotation_for_chr10_22"
$BEDOPS_PATH/bedmap --echo --echo-map-id --delim '\t' $partitioned_chr10_22_only $merged_annotation_chr10_22_only > $annotation_for_chr10_22

echo "combine them together!"
final_partitioned_file="final_partitioned_annotations.bed"
cat $worked_chromosome_partitioned_file > $final_partitioned_file
cat $annotation_for_chr10_22 >> $final_partitioned_file

# Finally, dump everything into MySQL database named: Gene_Annotations37/38
#
#printf "dump all partitioned exons into MySQL database Gene_Annotations37/38 ==> exonAnnotations.pl $final_partitioned_file \n"
$PERL_PATH $SCRIPT_PATH/exonAnnotations.pl "$final_partitioned_file" "$gene_db_version"

####
#END
####
