#!/bin/bash
# nodes=1:ppn=4,mem=32gb
SOURCE=`pwd`;
#JAVA="/hgsc_software/java/latest/bin/java";
JAVA="/hgsc_software/java/jdk1.8.0_74/bin/java";
IN_FILE="$1";
REDUCERATE="$2";
REFERENCE="$3";

## check the inputs
if [ "$#" -lt 2 ]; then
	echo "USAGE: downsample_bam_cram.sh   input_file   reduction_rate   reference_sequence(for cram file only)"
	exit;
fi;

##Check if the input file is present.
if [ ! -f ${IN_FILE} ]; then
    echo "Can't find $IN_FILE";
    exit;
fi;

##Check to see if the input file is a bam file
if [ ${IN_FILE: -4} == ".bam" ]; then
	#${JAVA} -Xmx15G -jar /hgsc_software/picard/picard-tools-1.54/DownsampleSam.jar I=${IN_FILE} O=${IN_FILE}_downsampled_${REDUCERATE}.bam P=${REDUCERATE};
	${JAVA} -Xmx32G -jar /hgsc_software/picard/picard-tools-2.6.0/picard.jar DownsampleSam I=${IN_FILE} O=${IN_FILE}_downsampled_${REDUCERATE}.bam P=${REDUCERATE} VALIDATION_STRINGENCY=SILENT TMP_DIR=/space1/tmp;

	${JAVA} -Xmx32G -jar /hgsc_software/picard/picard-tools-2.6.0/picard.jar MarkDuplicates TMP_DIR=/space1/tmp AS=TRUE M=/dev/null VALIDATION_STRINGENCY=SILENT I=${IN_FILE}_downsampled_${REDUCERATE}.bam O=${IN_FILE}_downsampled_${REDUCERATE}.marked.bam

	/hgsc_software/samtools/latest/samtools index ${IN_FILE}_downsampled_${REDUCERATE}.marked.bam
fi;

##check to see if the input file is a cram file
if [ ${IN_FILE: -5} == ".cram" ]; then
	## check if reference is available
	if [[ -z ${REFERENCE} ]]; then
		echo "Reference is not available for cram file. Exiting...";
		exit;
	fi;

	echo "processing cram file!"
	/hgsc_software/samtools/latest/samtools view -h -C -s ${REDUCERATE} ${IN_FILE} -o ${IN_FILE}_downsampled_${REDUCERATE}.cram;

	/hgsc_software/java/jdk1.8.0_74/bin/java -Xmx32G -jar /stornext/snfs5/next-gen/scratch/phuang/software/picard-2.10.7/build/libs/picard.jar MarkDuplicates TMP_DIR=/space1/tmp AS=TRUE M=/dev/null VALIDATION_STRINGENCY=SILENT I=${IN_FILE}_downsampled_${REDUCERATE}.cram O=${IN_FILE}_downsampled_${REDUCERATE}.marked.bam R=${REFERENCE}

	/hgsc_software/samtools/latest/samtools index ${IN_FILE}_downsampled_${REDUCERATE}.marked.bam
fi;
