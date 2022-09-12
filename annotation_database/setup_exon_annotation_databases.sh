#!/bin/bash

# This script is used to setup the gene/exon annotation database for sequencing analysis
# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 2 ]]; then
	echo "Illegal Number of Parameters"
	echo "/SCRIPT_PATH/setup_exon_annotation_databases.sh output_directory db_version(hg38 or hg37)"
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

# get the gene annotation version
gene_db_version=$2

# now process miRNA file
new_miRNA_file=hsa_miRNA_simplified.bed

if [ "$gene_db_version" == "hg38" ]; then
    wget https://mirbase.org/ftp/CURRENT/genomes/hsa.gff3
	awk -F "\t|" '{print $1"\t"$4"\t"$5"\t"$9}' hsa.gff3 | grep '^chr' | tr -s ';' '\t' | cut -f1,2,3,4,6 | sed s/ID=// | sed s/Name=// | awk '{t=$5; $5=$4; $4=t; print}' | awk '{ $4=$4"|exon_1|miRNA_1="$5; print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > $new_miRNA_file
	#printf "Producing simplified miRNA file $new_miRNA_file \n"
else 
    wget https://mirbase.org/ftp/20/genomes/hsa.gff3
	awk -F "\t| " '{print $1"\t"$4"\t"$5"\t"$9}' hsa.gff3 | grep '^chr' | tr -s ';' '\t' | cut -f1,2,3,4,6 | sed s/ID=// | sed s/Name=// | sed s/^chr//gi | awk '{t=$5"\t"; $5=$4"\t"; $4=t"\t"; print}' | awk '{ $4=$4"|exon_1|miRNA_1="$5; print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > $new_miRNA_file
	#new_miRNA_file=`basename ${miRNA_FILE}`_rearranged.bed
	#$SCRIPT_PATH/rearrange_miRNA_bed.py -i $miRNA_FILE > $new_miRNA_file
	#printf "Producing rearranged miRNA file $new_miRNA_file \n"
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

#if [ ${HGNC: -4} == ".txt" ]; then
#	mv $HGNC "hgnc_complete_set"
#fi

# now we need to download the newest gene annotation from various sources (such as refseq, ccds and gencode)
#
if [ "$gene_db_version" == "hg38" ]; then
	echo "For hg38"
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz .
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ccdsGene.txt.gz .
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV41.txt.gz .
	#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeCompV41.txt.gz .
else
	echo "For hg19"
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz .
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ccdsGene.txt.gz .
	rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV41lift37.txt.gz .
	#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeCompV41lift37.txt.gz  .
fi

# untar the zipped files
# it seems we need to use the CDS start/end for the final calculation, so we will need to add them here
#
if [ "$(uname)" == "Darwin" ]; then
	echo "For Darwin"
    ls $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
	if [ "$gene_db_version" == "hg38" ]; then
		ls $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; $SCRIPT_PATH/generate_bed_file.pl "$FILE.tmp" "$gene_db_version" > "$FILE.bed"; done;
	else
		ls $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; $SCRIPT_PATH/generate_bed_file.pl "$FILE.tmp" "$gene_db_version" | sed s/^chr//gi >> "$FILE.bed"; done;
	fi
elif [ "$(uname)" == "Linux" ]; then
	echo "For Linux"
    ls --color=never $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
	if [ "$gene_db_version" == "hg38" ]; then
		ls --color=never $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; $SCRIPT_PATH/generate_bed_file.pl "$FILE.tmp" "$gene_db_version" > "$FILE.bed"; done;
	else
		ls --color=never $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; $SCRIPT_PATH/generate_bed_file.pl "$FILE.tmp" "$gene_db_version" | sed s/^chr//gi > "$FILE.bed"; done;
	fi
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
	/stornext/snfs130/NGIRD/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$tmp_sorted_bed" | uniq - > $tmp_sorted_merged_bed

	# dump everything into a database
	# The -i option of grep says to ignore case.
	# The -q option says to not emit output and exit after the first match.
	# The -F option says to treat the argument as a string rather than a regular expression.
	#
	if echo $f | grep -iqF "ref" ; then
	#if [[ $f == *"ref"* ]]; then
		echo "dumping RefSeq exons"
		#$SCRIPT_PATH/dumpExonAnnotationFromIndividualSource.pl $tmp_sorted_merged_bed $gene_db_version "refseq"
	elif echo $f | grep -iqF "ccds"; then
	#elif [[ $f == *"ccds"* ]]; then
		echo "dumping CCDS exons"
		#$SCRIPT_PATH/dumpExonAnnotationFromIndividualSource.pl $tmp_sorted_merged_bed $gene_db_version "ccds"
	elif echo $f | grep -iqF "gen"; then
	#elif [[ $f == *"gen"* ]] ; then
		echo "dumping Gencode exons"
		#$SCRIPT_PATH/dumpExonAnnotationFromIndividualSource.pl $tmp_sorted_merged_bed $gene_db_version "gencode"
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
/stornext/snfs130/NGIRD/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$combined_sorted_bed_file" | uniq - > "$combined_sorted_bed_file_merged_perfect_matches"

# now need to extract all the exon info for all genes for batch analysis
#
#echo "running /hgsc_software/perl/perl-5.18.2/bin/perl $SCRIPT_PATH/geneExonAnnotationForBatchAnalysis.pl $combined_sorted_bed_file_merged_perfect_matches $gene_db_version"
#/hgsc_software/perl/perl-5.18.2/bin/perl $SCRIPT_PATH/geneExonAnnotationForBatchAnalysis.pl "$combined_sorted_bed_file_merged_perfect_matches" "$gene_db_version"

# Now it is the time to generate exon partition file
#
exons_intersected="$combined_sorted_bed_file"_intersected
exons_partitioned="$combined_sorted_bed_file"_merged_and_partitioned
printf "partition $combined_sorted_bed_file_merged_perfect_matches to produce $exons_partitioned using bedops \n"
#/stornext/snfs130/NGIRD/scratch/phuang/software/bin/bedops -p "$combined_sorted_bed_file_merged_perfect_matches" | /stornext/snfs130/NGIRD/scratch/phuang/software/bin/bedmap	--echo	--echo-map-id --delim '\t' - $combined_sorted_bed_file_merged_perfect_matches | uniq - > "$exons_partitioned"
/stornext/snfs130/NGIRD/scratch/phuang/software/bin/bedops -p "$combined_sorted_bed_file_merged_perfect_matches" | /hgsc_software/BEDTools/bedtools-2.26.0/bin/bedtools intersect -a - -b "$combined_sorted_bed_file_merged_perfect_matches"  -wb | cut -f1-3,7 > $exons_intersected
/hgsc_software/BEDTools/bedtools-2.26.0/bin/bedtools merge -i $exons_intersected -d -1 -c 4 -o collapse > "$exons_partitioned"
exons_partitioned_fixed="combined_sorted_bed_file"_merged_and_partitioned_fixed
sed 's/,/;/g' "$exons_partitioned" > "$exons_partitioned_fixed"

# Finally, dump everything into MySQL database named: Gene_Annotations37/38
#
printf "dump all partitioned exons into MySQL database Gene_Annotations37/38 ==> exonAnnotations.pl $final_partitioned_file \n"
/hgsc_software/perl/perl-5.18.2/bin/perl $SCRIPT_PATH/exonAnnotations.pl "$exons_partitioned_fixed" "$gene_db_version"

####
#END
####
