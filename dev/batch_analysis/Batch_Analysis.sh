#!/bin/bash

# this script is used to setup the Batch run for Batch analysis of sequencing data from the output of ExCID

## check the inputs
if [[ "$#" -lt 12  ||  "$#" -gt 14 ]]; then
	echo "USAGE: Batch_Analysis.sh [options]"
	echo -e "\t-a: The integer value from 1-5 to specify which annotation source you are interested in (Default 1)";
	echo -e "\t\t: 1: RefSeq annotation only (Default)";
	echo -e "\t\t: 2: CCDS annotation only";
	echo -e "\t\t: 3: Vega annotation only";
	echo -e "\t\t: 4: Gencode annotation only";
	echo -e "\t\t: 5: combined (all of the above)";
	echo -e "\t-b: the file name that contains gene exons to be examined in bed format (either VCRome+PKv2 or user defined) [Mandatory]";
   	echo -e "\t-d: the DB Version (such as hg19 or hg38) [Mandatory]";
	echo -e "\t-f: the input file that contains a list of low coverage reporting file names (absolute path) [Mandatory]";
	echo -e "\t-g: the file name (absolute path) that contains a list of genes user interested in";
   	echo -e "\t-o: the absolute path of your output directory [Mandatory]"
   	echo -e "\t-p: the Percentage of samples have low coverage (such as 30 or 50) to be reported [Mandatory]";
	exit;
fi;

# Use -gt 1 to consume two arguments per pass in the loop (eg: each argument has a corresponding value to go with it).
# Use -gt 0 to consume one or more arguments per pass in the loop (eg: some arguments don't have corresponding values to go with it)
# the --default option will handle such cases (ie, no corresponding values)
while [[ $# -gt 1 ]]; do
	key="$1"

	case $key in 
		-f|--file)
		# this file contains a list of low coverage report files in absolute path
		IN_FILE="$2"
		shift # pass argument
		;;
		-p|--percentage)
		Percentage="$2"
		shift # pass argument
		;;
		-a|--annotation)
		Annotation_Source="$2"
		shift
		;;
		-b|--bed)
		BED_FILE="$2"
		shift
		;;
		-d|--database)
		Database_Version="$2"
		shift # pass argument
		;;
		-g|--gene)
		GENE_LIST_FILE="$2"
		shift
		;;
		-o|--outputdir)
		OUTPUTDIR="$2"
		shift
		;;
		--default)
		DEFAULT=YES
		;;
		*)
		echo "you have entered an unknown option"       # for any unknown option
		;;
	esac
	shift # past argument or value
done

# check if Annotation_Source is set, if not, set it to '1'
# "-z" tests for a zero-length string.
# the inverse of -z is -n if [ -n "$VAR" ]; 
#
if [[ -z "$Annotation_Source" ]]; then
	echo "Annotation_Source is not set. Hence use default: 1 for RefSeq"
	Annotation_Source='1'
else
	echo "Annotation_Source set!"
fi

# calculate the threshold
FILE_TOTAl=`less $IN_FILE | wc -l`
THRESHOLD_f=$(($Percentage * $FILE_TOTAl))
THRESHOLD_f=$(($THRESHOLD_f/100))
THRESHOLD_i=`echo ${THRESHOLD_f%%.*}`
echo "threshold float is $THRESHOLD_f"
echo "threshold int is $THRESHOLD_i"

## open the input file list to see if they exist
for FILE in `cat $IN_FILE`; do
	if [ ! -f ${FILE} ]; then
		echo "Can't find $FILE";
		exit;
	fi;
done;

## need to convert input files into files without comments and concatenate all tmp file names together
SAMPLE_LIST=""
for FILE in `cat $IN_FILE`; do egrep -v "^#|^SeqStats" $FILE | sort -k1,1 -V -s -k2,2n -k3,3n > $OUTPUTDIR/$(basename $FILE).tmp; SAMPLE_LIST=$(printf "%s %s" "$SAMPLE_LIST" "$OUTPUTDIR/$(basename ${FILE}).tmp"); done
#echo "File list:"
#echo $SAMPLE_LIST

## intersect all of them (low coverage files) together using bedtools multiinter
#
echo "Now Intersect all the low coverage files"
LOW_COVERAGE_FILE="Batch_Low_COVERAGE_Multiinter.bed"
/hgsc_software/BEDTools/latest/bin/bedtools multiinter -i $SAMPLE_LIST | awk -F "\t" -v threshold="$THRESHOLD_i" ' $4>=threshold {print $1 "\t" $2 "\t" $3 "\t" $3-$2;} ' > $OUTPUTDIR/$LOW_COVERAGE_FILE 

## merge low coverage regions if they overlap
#
echo "Merge the low coverage regions"
LOW_COVERAGE_Merged_File="Batch_Low_COVERAGE_Merged.bed"
/hgsc_software/BEDTools/latest/bin/bedtools merge -i $OUTPUTDIR/$LOW_COVERAGE_FILE > $OUTPUTDIR/$LOW_COVERAGE_Merged_File

## now dump the official merged exon annotation (from MySQL database) into a bed file
#
echo "Fetch the official exon information"
/stornext/snfs5/next-gen/scratch/phuang/git_repo/dev/batch_analysis/src/batch_analysis -f $OUTPUTDIR/$LOW_COVERAGE_Merged_File -d $Database_Version -o $OUTPUTDIR -m 1 -a $Annotation_Source

## sort the official annotation bed file
#
echo "Now sort merged exons files"
sort_flag="-k1,1 -V -s -k2,2n"
Merged_Exons="Official_Merged_Exons.bed"
Merged_Exons_Sorted="Official_Merged_Exons_Sorted.bed"
sort $sort_flag $OUTPUTDIR/$Merged_Exons > $OUTPUTDIR/$Merged_Exons_Sorted

## intersect merged_exon annotations with $LOW_COVERAGE_Merged_File
#
echo "Intersect the official exon annotation with the low coverage regions"
Low_Coverage_Merged_Annotation=$LOW_COVERAGE_Merged_File"_Annotations"
/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $OUTPUTDIR/$Merged_Exons_Sorted -b $OUTPUTDIR/$LOW_COVERAGE_Merged_File -wb > $OUTPUTDIR/$Low_Coverage_Merged_Annotation

## now intersect the merged_exon annotations with the input bed file to check which exons are NOT covered!
#
echo "Find out which official exon annotations are not covered by the Capture design using intersect -v option"
NOT_Covered_Exons="NOT_Covered_Exons.bed"
/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $OUTPUTDIR/$Merged_Exons_Sorted -b $BED_FILE -v > $OUTPUTDIR/$NOT_Covered_Exons

## now need to calculate the detailed information regarding gene annotation within low coverage regions
## read in our design not covered exons (pass in as an argument)
##
echo "Obtain the annotation"
/stornext/snfs5/next-gen/scratch/phuang/git_repo/dev/batch_analysis/src/batch_analysis -f $OUTPUTDIR/$Low_Coverage_Merged_Annotation -d $Database_Version -o $OUTPUTDIR -g $GENE_LIST_FILE -m 2 -u $OUTPUTDIR/$NOT_Covered_Exons -a $Annotation_Source
