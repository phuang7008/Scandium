#!/bin/bash

# This script is used to setup the gene cds coordinate database for the percentage analysis
# Here we are only interested in RefSeq CDS coordinates. 
# For CCDS or Gencode, we will have to modify code accordingly
# YES, Qiaoyan is interested in Gencode annotation as well
# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 3 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/setup_gene_coords_databases.sh output_dir db_version(hg38 or hg37) db_name"
	exit
fi

# get the working directory and cd to the directory
#
BASEDIR=$1
cd $BASEDIR
pwd
#echo `$BASEDIR`
printf "$BASEDIR\n"

gene_db_version=$2
db_name=$3

# use case insentive comparison
shopt -s nocasematch

renamed_file=""
project_name=""

if [[ $db_name == "refseq" ]]; then
	# here we are only interested in refseq source of annotation/coords.
	#
	if [ "$gene_db_version" == "hg38" ]; then
		rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz .
	else
		rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz .
	fi

	renamed_file='RefSeq_renamed'
	project_name="RefSeq"

else
	# Gencode
	#
	if [ "$gene_db_version" == "hg38" ]; then
		rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV28.txt.gz .
		#rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeCompV28.txt.gz .
	else
		rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeBasicV19.txt.gz .
	fi

	renamed_file='Gencode_renamed'
	project_name="Gencode"
fi

# untar the zipped files
# several refseq sequences could have the same name, even if they refer to different genomic locations
# For example: The name will be changed from NM_000015 to NM_000015-1, NM_000015-2 etc
#
if [ "$(uname)" == "Darwin" ]; then
    ls $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/renameRefSeq.pl "$FILE.tmp" > $renamed_file; done;
elif [ "$(uname)" == "Linux" ]; then
    ls --color=never $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls --color=never $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/renameRefSeq.pl "$FILE.tmp" > $renamed_file; done;
fi

bed_file="$renamed_file"_bed
/stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/generate_refseq_exon_bed.py -i "$renamed_file" -v $gene_db_version > "$bed_file"

# save all the .txt.gz downloaded source files in case we need it in the future!
#
tmp_timestamp=$(date +%Y'_'%m'_'%d)
printf "time stamp is $tmp_timestamp\n"
for f in $BASEDIR/*.txt; do gzip $f; mv "$f.gz" "$f.$tmp_timestamp.gz"; done

# Before we continue, we need to sort the bed_files 
#
sorted_bed_file="$bed_file"_sorted
echo "sort $bed_file to produce $sorted_bed_file"
sort_flag="-k1,1 -V -s -k2,2n -k3,3n"
sort $sort_flag -o $sorted_bed_file $bed_file

# Next, we need to merge those exons with perfect coordinates.
#
sorted_bed_file_merged_perfect_matches="$sorted_bed_file"_merged_perfect_matches
echo "merge perfect matched regions in $sorted_bed_file to produce $sorted_bed_file_merged_perfect_matches"
/stornext/snfs5/next-gen/scratch/phuang/software/bin/bedmap --echo-map-range --echo-map-id --delim '\t' --fraction-both 1.0 "$sorted_bed_file" | uniq - > "$sorted_bed_file_merged_perfect_matches"

# to better scale the analysis for other project such as eMerge and CCDS, we will dump CDS here and do the rest when we are running it
#
#echo "dump all the CDS into the database"
#/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/refseqCDSAnnotation.pl "$sorted_bed_file_merged_perfect_matches" "$gene_db_version"

#echo "DONE for CDS into the database!"

# dump everything into MySQL table named: Gene_DB-Name_Coords37/38 (either RefSeq or Gencode or some other databases)
#
echo "dump all exons into MySQL database Gene_RefSeq_Coords37/38 ==> fetchGeneCDSCoords.pl $exon_target_intersect_file"
/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/annotation_database/fetchGeneCDSCoords.pl "$sorted_bed_file_merged_perfect_matches" "$gene_db_version" "$project_name"

shopt -u nocasematch

####
#END
####
