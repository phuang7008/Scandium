#!/bin/bash

# this script is used to setup the Batch run for Batch analysis of sequencing data from the output of ExCID

## check the inputs
function help() {
	echo "USAGE: Batch_Analysis.sh [options]"
	echo -e "\t-a: The integer value from 1-5 to specify which annotation source you are interested in (Default 1)";
	echo -e "\t\t: 1: RefSeq annotation only (Default)";
	echo -e "\t\t: 2: CCDS annotation only";
	echo -e "\t\t: 3: Vega annotation only";
	echo -e "\t\t: 4: Gencode annotation only";
	echo -e "\t\t: 5: combined (all of the above)";
	echo -e "\t-b: the target bed file that contains exons to be examined (either VCRome+PKv2 or user defined) [Mandatory]";
   	echo -e "\t-d: the DB Version (such as hg19 or hg38) (Default hg38)";
	echo -e "\t-i: the input file that contains a list of low coverage reporting file names (absolute path) [Mandatory]";
	echo -e "\t-f: the file that contains user defined database";
	echo -e "\t-g: the file name (absolute path) that contains a list of genes user interested in";
   	echo -e "\t-o: the absolute path of your output directory [Mandatory]"
   	echo -e "\t-p: the Percentage of samples have low coverage (such as 10 or 35) [Mandatory]";
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
if [[ (-z "$OUTPUTDIR" ) || (-z "$IN_FILE") || (-z "$Percentage") || (-z "$BED_FILE") ]]; then
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
for FILE in `cat $IN_FILE`; do grep -v "^#" $FILE | cut -f1,2,3 | /hgsc_software/BEDTools/latest/bin/bedtools sort -i - > $OUTPUTDIR/$(basename $FILE).tmp; SAMPLE_LIST=$(printf "%s %s" "$SAMPLE_LIST" "$OUTPUTDIR/$(basename ${FILE}).tmp"); done
#echo "File list:"
#echo $SAMPLE_LIST
echo

## intersect all of them (low coverage files) together using bedtools multiinter
echo "Intersect all sample files using bedtools multiinter"
LOW_COVERAGE_FILE="Batch_Low_COVERAGE_Multiinter.bed"
/hgsc_software/BEDTools/latest/bin/bedtools multiinter -i $SAMPLE_LIST | awk -F "\t" -v threshold="$THRESHOLD_i" ' $4>=threshold {print $1 "\t" $2 "\t" $3 "\t" $3-$2;} ' > $OUTPUTDIR/$LOW_COVERAGE_FILE

## merge low coverage regions if they overlap
echo "Merged low coverage regions if they overlap"
LOW_COVERAGE_Merged_File="Batch_Low_COVERAGE_Merged.bed"
/hgsc_software/BEDTools/latest/bin/bedtools sort -i $OUTPUTDIR/$LOW_COVERAGE_FILE | /hgsc_software/BEDTools/latest/bin/bedtools merge -i - > $OUTPUTDIR/$LOW_COVERAGE_Merged_File

# After talking to Qiaoyan, we decided not to use MySQL database anymore ...
# But I think it might be helpful in the future to use it! So I added it in
## now dump the official exon annotation (from MySQL database: default RefSeq) into a bed file
if [ -z "$Annotation_Source" ]; then
	Annotation_Source=1
fi

## intersect exon annotations with $LOW_COVERAGE_Merged_File
Low_Coverage_Merged_Annotation=$LOW_COVERAGE_Merged_File"_Annotations"

if [ -z "$User_Defined_DB" ]; then
	# no user-defined database (annotations), so let's get the Official exon annotation from MySQL database
	#
	echo "Dump the official exon annotation from MySQL database) into a bed file"
	/stornext/snfs5/next-gen/scratch/phuang/dev/batch_analysis/src/batch_analysis -i $OUTPUTDIR/$LOW_COVERAGE_Merged_File -d $Database_Version -o $OUTPUTDIR -m 1 -a $Annotation_Source -t $BED_FILE

	## sort the official annotation bed file
	#echo "Sort the official annotation bed file"
	Official_Exons="Official_Exons.bed"
	#sort -k1,1 -V -s $Official_Exons > $Merged_Exons_Sorted

	# MySQL Official database is too big, we need to remove those un-covered exons
	echo "Intersect official exon annotation with target bed file to remove un-covered exons."
	Targeted_Exon_Annotation="Targeted_Exon_Annotation"
	/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $BED_FILE -b $Official_Exons -wao | awk -F"\t" '($7!="") {print}' | cut -f1,2,3,7 > $Targeted_Exon_Annotation

	# Intersect the shortened Official database (ie, Targeted_Exon_Annotation)  with $LOW_COVERAGE_Merged_File
	/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $OUTPUTDIR/$LOW_COVERAGE_Merged_File -b $OUTPUTDIR/$Targeted_Exon_Annotation -wao > $Low_Coverage_Merged_Annotation

	echo "Now run the batch analysis"
	echo "/stornext/snfs5/next-gen/scratch/phuang/dev/batch_analysis/src/batch_analysis -i $OUTPUTDIR/$Low_Coverage_Merged_Annotation -d $Database_Version -o $OUTPUTDIR $GENE_LIST_FILE -m 2 -a $Annotation_Source -f $Targeted_Exon_Annotation -t $BED_FILE"

	/stornext/snfs5/next-gen/scratch/phuang/dev/batch_analysis/src/batch_analysis -i $OUTPUTDIR/$Low_Coverage_Merged_Annotation -d $Database_Version -o $OUTPUTDIR $GENE_LIST_FILE -m 2 -a $Annotation_Source -f $Targeted_Exon_Annotation -t $BED_FILE

else

	# Intersect the user-defined database/annotation with $LOW_COVERAGE_Merged_File
	echo "Intersect user-defined annotations with $LOW_COVERAGE_Merged_File"
	/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $OUTPUTDIR/$LOW_COVERAGE_Merged_File -b $User_Defined_DB -wao > $Low_Coverage_Merged_Annotation

	echo "Now run the batch analysis"
	echo "/stornext/snfs5/next-gen/scratch/phuang/dev/batch_analysis/src/batch_analysis -i $OUTPUTDIR/$Low_Coverage_Merged_Annotation -d $Database_Version -o $OUTPUTDIR $GENE_LIST_FILE -m 2 -a $Annotation_Source -f $User_Defined_DB -t $BED_FILE"

	/stornext/snfs5/next-gen/scratch/phuang/dev/batch_analysis/src/batch_analysis -i $OUTPUTDIR/$Low_Coverage_Merged_Annotation -d $Database_Version -o $OUTPUTDIR $GENE_LIST_FILE -m 2 -a $Annotation_Source -f $User_Defined_DB -t $BED_FILE
fi 

# now remove all the tmp files we just created
echo "Remove all tmp files just created!"
for FILE in `cat $IN_FILE`; do rm -f $OUTPUTDIR/$(basename $FILE).tmp; done

exit
