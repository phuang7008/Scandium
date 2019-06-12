#!/bin/bash
# nodes=1:ppn=4,mem=36gb
SOURCE=`pwd`;
#JAVA="/hgsc_software/java/latest/bin/java";
JAVA="/hgsc_software/java/jdk1.8.0_74/bin/java";

## check the inputs
if [ "$#" -lt 3 ]; then
	echo "USAGE: downsample_bam_cram.sh   input_file   reduction_rate   output_dir  reference_sequence(for cram file only)"
	exit;
fi;

IN_FILE=$1;
REDUCERATE=$2;
output_dir=$3;
REFERENCE=$4;

base_filename=`basename $IN_FILE`
downsampled_filename=${output_dir}/${base_filename}_downsampled_${REDUCERATE}.bam
marked_bam=${output_dir}/${base_filename}_downsampled_${REDUCERATE}.marked.bam

##Check if the input file is present.
if [ ! -f ${IN_FILE} ]; then
    echo "Can't find $IN_FILE";
    exit;
fi;

##Check to see if the input file is a bam file
if [ ${IN_FILE: -4} == ".bam" ]; then
	${JAVA} -jar /hgsc_software/picard/picard-tools-2.20.1/picard.jar DownsampleSam I=${IN_FILE} O=$downsampled_filename P=${REDUCERATE} VALIDATION_STRINGENCY=SILENT TMP_DIR=/space1/tmp
	#${JAVA} -Xmx36G -jar /hgsc_software/picard/picard-tools-2.20.1/picard.jar DownsampleSam I=${IN_FILE} O=${downsampled_filename} P=${REDUCERATE} VALIDATION_STRINGENCY=SILENT TMP_DIR=/space1/tmp

	${JAVA} -jar /hgsc_software/picard/picard-tools-2.20.1/picard.jar MarkDuplicates TMP_DIR=/space1/tmp AS=TRUE M=/dev/null VALIDATION_STRINGENCY=SILENT I=$downsampled_filename O=$marked_bam
	#${JAVA} -Xmx36G -jar /hgsc_software/picard/picard-tools-2.20.1/picard.jar MarkDuplicates TMP_DIR=/space1/tmp AS=TRUE M=/dev/null VALIDATION_STRINGENCY=SILENT I=${downsampled_filename} O=${marked_bam}

	/hgsc_software/samtools/latest/samtools index $marked_bam
fi;

##check to see if the input file is a cram file
if [ ${IN_FILE: -5} == ".cram" ]; then
	## check if reference is available
	if [[ -z ${REFERENCE} ]]; then
		echo "Reference is not available for cram file. Exiting...";
		exit;
	fi;

	echo "processing cram file!"
    ${JAVA} -jar /hgsc_software/picard/picard-tools-2.20.1/picard.jar DownsampleSam I=${IN_FILE} O=$downsampled_filename P=${REDUCERATE} VALIDATION_STRINGENCY=SILENT TMP_DIR=/space1/tmp REFERENCE_SEQUENCE=${REFERENCE}
    #${JAVA} -Xmx36G -jar /hgsc_software/picard/picard-tools-2.20.1/picard.jar DownsampleSam I=${IN_FILE} O=$downsampled_filename P=${REDUCERATE} VALIDATION_STRINGENCY=SILENT TMP_DIR=/space1/tmp REFERENCE_SEQUENCE=${REFERENCE}

	${JAVA} -jar /hgsc_software/picard/picard-tools-2.20.1/picard.jar MarkDuplicates TMP_DIR=/space1/tmp AS=TRUE M=/dev/null VALIDATION_STRINGENCY=SILENT I=$downsampled_filename O=$marked_bam REFERENCE_SEQUENCE=${REFERENCE}
	#${JAVA} -Xmx36G -jar /hgsc_software/picard/picard-tools-2.20.1/picard.jar MarkDuplicates TMP_DIR=/space1/tmp AS=TRUE M=/dev/null VALIDATION_STRINGENCY=SILENT I=${downsampled_filename} O=${marked_bam} REFERENCE_SEQUENCE=${REFERENCE}

	/hgsc_software/samtools/latest/samtools index ${marked_bam}
fi;
