#!/bin/bash

# This script is used to setup the gene/exon annotation database for sequencing analysis
# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 3 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/setup_exon_annotation_databases.sh miRNA_file output_directory db_version(hg38 or hg37)"
	exit
fi

# get the working directory and cd to the directory
BASEDIR=$2
cd $BASEDIR
pwd
#echo `$BASEDIR`
printf "$BASEDIR\n"

# get the gene annotation version
gene_db_version=$3

# since miRNA file can't be downloaded using rsync or ftp, we have to download it first and tell the script where to find it
# things might change in the future and we might be able to use rsync or ftp in the future! 
# But for now, we need to do it manually!
miRNA_FILE=$1

# now process miRNA file
new_miRNA_file=""

if [ "$gene_db_version" == "hg38" ]; then
	new_miRNA_file=`basename ${miRNA_FILE}`_simplified.bed
	awk -F "\t| " '{print $1"\t"$4"\t"$5"\t"$9}' $miRNA_FILE | grep 'chr' | tr -s ';' '\t' | cut -f1,2,3,4,6 | sed s/ID=// | sed s/Name=// | awk '{t=$5; $5=$4; $4=t; print}' | awk '{ $4=$4"\|exon_1\|miRNA_1="$5; print}' > $new_miRNA_file
	#awk -F "\t| " '{print $1"\t"$4"\t"$5"\t"$9}' $miRNA_FILE | grep 'chr' | tr -s ';' '\t' | cut -f1,2,3,4,6 | sed s/ID=// | sed s/Name=// | sed s/chr//gi | awk '{t=$5; $5=$4; $4=t; print}' | awk '{ $4=$4"_exon_0_1="$5; print}' > $new_miRNA_file
	printf "Producing simplified miRNA file $new_miRNA_file \n"
else 
	new_miRNA_file=`basename ${miRNA_FILE}`_rearranged.bed
	/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/rearrange_miRNA_bed.py -i $miRNA_FILE > $new_miRNA_file
	printf "Producing rearranged miRNA file $new_miRNA_file \n"
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
#/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/processHGNCtoDB.pl "$HGNC_simplified" "$gene_db_version"

# now get rid of .txt file extension before preceed, as we will unzip the tar file into .txt
#
for f in $BASEDIR/*.txt; do mv -- "$f" "${f%.txt}"; done

#if [ ${HGNC: -4} == ".txt" ]; then
#	mv $HGNC "hgnc_complete_set"
#fi

# now we need to download the newest gene annotation from various sources (such as refseq, ccds and gencode)
#
if [ "$gene_db_version" == "hg38" ]; then
	echo "For hg38"
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz .
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ccdsGene.txt.gz .
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV26.txt.gz .
	#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeCompV26.txt.gz .
	#rsync -a -P rsync://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3 .
else
	echo "For hg19"
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz .
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ccdsGene.txt.gz .
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vegaGene.txt.gz .
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz .
	#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeCompV19.txt.gz .
	#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vegaGtp.txt.gz .
	#rsync -a -P rsync://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3 .
fi

# untar the zipped files
# it seems we need to use the CDS start/end for the final calculation, so we will need to add them here
#
if [ "$(uname)" == "Darwin" ]; then
	echo "For Darwin"
    ls $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/generate_bed_file.pl "$FILE.tmp" "$gene_db_version" > "$FILE.bed"; done;
elif [ "$(uname)" == "Linux" ]; then
	echo "For Linux"
    ls --color=never $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls --color=never $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/generate_bed_file.pl "$FILE.tmp" "$gene_db_version" > "$FILE.bed"; done;
fi

# remove all the tmp files files
#
for f in $BASEDIR/*.tmp; do rm -f $f; done
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
	/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$tmp_sorted_bed" | uniq - > $tmp_sorted_merged_bed

	# dump everything into a database
	# The -i option of grep says to ignore case.
	# The -q option says to not emit output and exit after the first match.
	# The -F option says to treat the argument as a string rather than a regular expression.
	#
	if echo $f | grep -iqF "ref" ; then
	#if [[ $f == *"ref"* ]]; then
		echo "dumping RefSeq exons"
		#/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/dumpExonAnnotationFromIndividualSource.pl $tmp_sorted_merged_bed $gene_db_version "refseq"
	elif echo $f | grep -iqF "ccds"; then
	#elif [[ $f == *"ccds"* ]]; then
		echo "dumping CCDS exons"
		#/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/dumpExonAnnotationFromIndividualSource.pl $tmp_sorted_merged_bed $gene_db_version "ccds"
	elif echo $f | grep -iqF "vega"; then
	#elif [[ $f == *"vega"* ]]; then
		echo "dumping VEGA exons"
		#/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/dumpExonAnnotationFromIndividualSource.pl $tmp_sorted_merged_bed $gene_db_version "vega"
	elif echo $f | grep -iqF "gen"; then
	#elif [[ $f == *"gen"* ]] ; then
		echo "dumping Gencode exons"
		#/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/dumpExonAnnotationFromIndividualSource.pl $tmp_sorted_merged_bed $gene_db_version "gencode"
	fi
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
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$combined_sorted_bed_file" | uniq - > "$combined_sorted_bed_file_merged_perfect_matches"

# now need to extract all the exon info for all genes for batch analysis
#
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/geneExonAnnotationForBatchAnalysis.pl "$combined_sorted_bed_file_merged_perfect_matches" "$gene_db_version"

# Now it is the time to generate exon partition file
#
exons_partitioned="$combined_sorted_bed_file"_merged_and_partitioned
printf "partition $combined_sorted_bed_file_merged_perfect_matches to produce $exons_partitioned using bedops \n"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedops -p "$combined_sorted_bed_file_merged_perfect_matches" | /stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap	--echo	--echo-map-id --delim '\t' - $combined_sorted_bed_file_merged_perfect_matches | uniq - > "$exons_partitioned"

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
awk -F"\t" '$1>9' $combined_sorted_bed_file_merged_perfect_matches | awk -F"\t" '$1<=22' > $merged_annotation_chr10_22_only

echo "now do the bedmap to generated annotated info for chr10 to chr22"
annotation_for_chr10_22="annotation_for_chr10_22"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo --echo-map-id --delim '\t' $partitioned_chr10_22_only $merged_annotation_chr10_22_only > $annotation_for_chr10_22

echo "combine them together!"
final_partitioned_file="final_partitioned_annotations.bed"
cat $worked_chromosome_partitioned_file > $final_partitioned_file
cat $annotation_for_chr10_22 >> $final_partitioned_file

# Finally, dump everything into MySQL database named: Gene_Annotations37/38
#
printf "dump all partitioned exons into MySQL database Gene_Annotations37/38 ==> exonAnnotations.pl $final_partitioned_file \n"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/exonAnnotations.pl "$final_partitioned_file" "$gene_db_version"

####
#END
####
