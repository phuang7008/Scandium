#!/bin/bash

# this script is used to setup the Batch run for Batch analysis of sequencing data from the output of ExCID

## check the inputs
function help() {
	echo "USAGE: Batch_Analysis.sh [options]"
	echo -e "\t-b: the target bed file that contains exon targets to be examined [Mandatory]";
   	echo -e "\t-d: the DB Version (such as hg19 or hg38) (Default hg38)";
	echo -e "\t-i: the input file that contains a list of low coverage reporting file names (absolute path) [Mandatory]";
	echo -e "\t-f: the file that contains user defined database";
	echo -e "\t-g: the file name (absolute path) that contains a list of genes user interested in";
   	echo -e "\t-o: the absolute path of your output directory [Mandatory]"
   	echo -e "\t-p: the Percentage of samples have low coverage (such as 10 or 35) [Mandatory]";
   	echo -e "\t-A: the batch_analysis program path [Mandatory]";
   	echo -e "\t-B: the bedtools program path [Mandatory]";
	exit;
}

# Use -gt 1 to consume two arguments per pass in the loop (eg: each argument has a corresponding value to go with it).
# Use -gt 0 to consume one or more arguments per pass in the loop (eg: some arguments don't have corresponding values to go with it)
# the --default option will handle such cases (ie, no corresponding values)
while [[ $# -gt 1 ]]; do
	key="$1"

	case $key in 
		-i|--inFile)
		# this file contains a list of low coverage report files in absolute path
		IN_FILE="$2"
		shift # pass argument
		;;
		-p|--percentage)
		Percentage="$2"
		shift # pass argument
		;;
		-b|--bed)
		BED_FILE="$2"
		shift
		;;
		-d|--database)
		Database_Version="$2"
		shift # pass argument
		;;
		-f|--udFile)
		User_Defined_DB="$2"
		shift
		;;
		-g|--gene)
		GENE_LIST_FILE="$2"
		shift
		;;
		-o|--outputdir)
		OUTPUTDIR="$2"
		shift
		;;
        -A|--batchanalysispath)
        BATCHANALYSIS="$2"
        shift
        ;;
        -B|--bedtoolspath)
        BEDTOOLS="$2"
        shift
        ;;
		--default)
		DEFAULT=YES
		;;
		*)
		echo $1
		echo $2
		echo "you have entered an unknown option"       # for any unknown option
		;;
	esac
	shift # past argument or value
done

# check user input first
#
if [[ (-z "$OUTPUTDIR" ) || (-z "$IN_FILE") || (-z "$Percentage") || (-z "$BED_FILE") || (-z "$BEDTOOLS") || (-z "$BATCHANALYSIS") ]]; then
	help
fi

if [ -z "$Database_Version" ]; then
	Database_Version="hg38"
fi

# handle gene list
if [ -z "$GENE_LIST_FILE" ]; then
	GENE_LIST_FILE=""
else
	GENE_LIST_FILE="-g $GENE_LIST_FILE"
fi

# calculate the threshold
FILE_TOTAl=`less $IN_FILE | wc -l`
THRESHOLD_f=$(($Percentage * $FILE_TOTAl))
THRESHOLD_f=$(awk "BEGIN {printf \"%0.2f\", ${THRESHOLD_f}/100}")
echo "threshold float is $THRESHOLD_f"

# usually we add 5 for ceiling, but this time, we will add 4 for ceiling
THRESHOLD_i=`perl -e "print int ($THRESHOLD_f + 0.4)"`
echo "threshold int is $THRESHOLD_i"

## open the input file list to see if they exist
for FILE in `cat $IN_FILE`; do
	if [ ! -f ${FILE} ]; then
		echo "Can't find $FILE";
		exit;
	fi;
done;

## need to convert input files into files without comments and concatenate all tmp file names together
## we also need to sort it as well
##
SAMPLE_LIST=""
for FILE in `cat $IN_FILE`; do grep -v "^#" $FILE | cut -f1,2,3 | $BEDTOOLS sort -i - > $OUTPUTDIR/$(basename $FILE).tmp; SAMPLE_LIST=$(printf "%s %s" "$SAMPLE_LIST" "$OUTPUTDIR/$(basename ${FILE}).tmp"); done
#echo "File list:"
#echo $SAMPLE_LIST
echo

## intersect all of them (low coverage files) together using bedtools multiinter
echo "Intersect all sample files using bedtools multiinter"
LOW_COVERAGE_FILE="Batch_Low_COVERAGE_Multiinter.bed"
$BEDTOOLS multiinter -i $SAMPLE_LIST | awk -F "\t" -v threshold="$THRESHOLD_i" ' $4>=threshold {print $1 "\t" $2 "\t" $3 "\t" $3-$2;} ' > $OUTPUTDIR/$LOW_COVERAGE_FILE

## merge low coverage regions if they overlap
echo "Merged low coverage regions if they overlap"
LOW_COVERAGE_Merged_File="Batch_Low_COVERAGE_Merged.bed"
$BEDTOOLS sort -i $OUTPUTDIR/$LOW_COVERAGE_FILE | $BEDTOOLS merge -i - > $OUTPUTDIR/$LOW_COVERAGE_Merged_File

# After talking to Qiaoyan, we decided not to use MySQL database anymore ...
#
if [ -z "$User_Defined_DB" ]; then
	# no user-defined database (annotations)
	#
	echo "$BATCHANALYSIS -i $OUTPUTDIR/$Low_Coverage_Merged_Annotation -d $Database_Version -o $OUTPUTDIR $GENE_LIST_FILE -t $BED_FILE"

    $BATCHANALYSIS -i $OUTPUTDIR/$Low_Coverage_Merged_Annotation -d $Database_Version -o $OUTPUTDIR $GENE_LIST_FILE -t $BED_FILE

else
	# Intersect with the user-defined database/annotation with $LOW_COVERAGE_Merged_File
    #
    Low_Coverage_Merged_Annotation=$LOW_COVERAGE_Merged_File"_Annotations"
	echo "Intersect user-defined annotations with $LOW_COVERAGE_Merged_File"
    $BEDTOOLS intersect -a $OUTPUTDIR/$LOW_COVERAGE_Merged_File -b $User_Defined_DB -wao > $OUTPUTDIR/$Low_Coverage_Merged_Annotation

	echo "Now run the batch analysis"
	echo "$BATCHANALYSIS -i $OUTPUTDIR/$Low_Coverage_Merged_Annotation -d $Database_Version -o $OUTPUTDIR $GENE_LIST_FILE -f $User_Defined_DB -t $BED_FILE"

    $BATCHANALYSIS -i $OUTPUTDIR/$Low_Coverage_Merged_Annotation -d $Database_Version -o $OUTPUTDIR $GENE_LIST_FILE -f $User_Defined_DB -t $BED_FILE
fi 

# now remove all the tmp files we just created
echo "Remove all tmp files just created!"
for FILE in `cat $IN_FILE`; do rm -f $OUTPUTDIR/$(basename $FILE).tmp; done

exit
