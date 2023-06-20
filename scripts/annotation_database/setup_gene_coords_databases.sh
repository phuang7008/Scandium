#!/bin/bash

# This script is used to setup the gene cds coordinate database for the percentage analysis
# Here we are only interested in 'refseq' or 'gencode' CDS annotations based on the user's option for db_name. 
# For CCDS or microRNA or any other database, users need to modify script accordingly

# First, we need to make sure that user has entered all the required parameters (or arguments)
#
if [[ $# -ne 6 ]]; then
	echo "Illegal Number of Parameters"
	echo "/SCRIPT_PATH/setup_gene_coords_databases.sh output_dir db_version(hg38 or hg37) db_name db_url python_interpreter perl_interpreter"
    echo "For example the gencode db_url: hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV43.txt.gz"
    echo "For example the refseq  db_url: hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz"
	exit
fi

# Don't want to hard code the path to the script. 
# Use the following to get the bash script path
#
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# get the working directory and cd to the directory
#
BASEDIR=$1
cd $BASEDIR
pwd
#echo `$BASEDIR`
printf "$BASEDIR\n"

gene_db_version=$2
db_name=$3
db_url=$4
python_interpreter=$5
perl_interpreter=$6

# use case insentive comparison
shopt -s nocasematch

renamed_file=""
project_name=""

rsync -a -P rsync://$db_url .

if [[ $db_name == "refseq" ]]; then
	# here we are only interested in refseq source of annotation/coords.
	#
	renamed_file='RefSeq_renamed'
	project_name="RefSeq"

elif [[ $db_name == "gencode" ]]; then
	# Gencode
	#
	renamed_file='Gencode_renamed'
	project_name="Gencode"
else
    echo "Please specify either 'refseq' or 'gencode' for processing"
    exit
fi

# untar the zipped files
# several refseq transcripts could have the same name, even if they refer to different genomic locations
# In this case, the name will be changed from NM_000015 to NM_000015-1, NM_000015-2 etc. using script renameRefSeq.pl
# but for gencode, this won't happen, the script just organize the file for the next step.
#
if [ "$(uname)" == "Darwin" ]; then
    ls $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; $perl_interpreter $SCRIPT_PATH/renameRefSeq.pl "$FILE.tmp" > $renamed_file; done;
elif [ "$(uname)" == "Linux" ]; then
    ls --color=never $BASEDIR/*.gz | while read FILE ; do gzip -d "$FILE" ; done ;
    ls --color=never $BASEDIR/*.txt | while read FILE ; do awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13}' "$FILE" > "$FILE.tmp"; $perl_interpreter $SCRIPT_PATH/renameRefSeq.pl "$FILE.tmp" > $renamed_file; done;
fi

bed_file="$renamed_file"_bed
$python_interpreter $SCRIPT_PATH/generate_annotation_bed.py -i "$BASEDIR/$renamed_file" -t $db_name -v $gene_db_version  > "$BASEDIR/$bed_file"

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

# dump everything into MySQL table named: Gene_DB-Name_Coords37/38 (either RefSeq or Gencode or some other databases)
#
echo "dump all exons into MySQL database Gene_RefSeq_Coords37/38 ==> fetchGeneCDSCoords.pl $exon_target_intersect_file"
#$perl_interpreter $SCRIPT_PATH/fetchGeneCDSCoords.pl "$sorted_bed_file" "$gene_db_version" "$project_name"

shopt -u nocasematch

####
#END
####
