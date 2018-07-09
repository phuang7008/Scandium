#!/bin/bash

# This script is used to facilitate the ExCID run
#

# First the input parameters
#
if [[ "$#" -lt 4 || "$#" -gt 31 ]]; then
	echo ""
	#/stornext/snfs5/next-gen/scratch/phuang/git_repo/dev/sequenceStats/coverage -h
	/stornext/snfs5/next-gen/scratch/phuang/dev/sequenceStats/coverage -h
	echo ""
	exit
fi;

# use -gt 0 to comsume 2 arguments by shifting twice for each case
# the --default option will handle such cases (ie, no corresponding values)
#
while [[ $# -gt 0 ]]; do
	key="$1"

	case $key in
		-i)
		INPUT_bam_FILE="$2"		# cram files are also accepted!
		shift # pass argument
		shift # pass value
		;;
		-o)
		OUTPUTDIR="$2"
		shift
		shift
		;;
		-b)
		MIN_BASE_QUAL="$2"
		shift
		shift
		;;
		-m)
		MIN_MAP_QUAL="$2"
		shift
		shift
		;;
		-n)
		Ns_REGION_FILE="$2"
		shift
		shift
		;;
		-p)
		PERCENTAGE="$2"
		shift
		shift
		;;
		-t)
		TARGET_FILE="$2"
		shift
		shift
		;;
		-D)
		DB_VERSION="$2"
		shift
		shift
		;;
		-H)
		HIGH_COVERAGE="$2"
		shift
		shift
		;;
		-L)
		LOW_COVERAGE="$2"
		shift
		shift
		;;
		-T)
		NUM_OF_THREADS="$2"
		shift
		shift
		;;
		-U)
		UPPER_BOUND_COVERAGE="$2"
		shift
		shift
		;;
		-a)
		ANNOTATION_ON='YES'
		shift
		;;
		-d)
		REMOVE_DUPLICATE_READS='YES'
		shift
		;;
		-s)
		REMOVE_SUPPLEMENTARY_READS='YES'
		shift
		;;
		-w)
		WRITE_WGS_COVERAGE='YES'
		shift
		;;
		-G)
		WRITE_WIG_FILE='YES'
		shift
		;;
		-W)
		WRITE_COV_FASTA_FILE='YES'
		shift
		;;
		--default)
		DEFAULT=YES
		;;
		*)
		echo "you have entered an unknown option"       # for any unknown option
		;;
	esac
done

# next, we will dump the Gene_RefSeq_CDS data out from the MySQL 
#
RefSeq_CDS_FILE="All_RefSeq_CDS"
/stornext/snfs5/next-gen/scratch/phuang/dev/sequenceStats/coverage -f $RefSeq_CDS_FILE

# do the intersect between RefSeq_CDS_FILE and input target bed file
#
exon_target_intersect_file="exon_target_intersect_for_gene_percentage_annotation"                            
echo "bedtools intersect to produce $exon_target_intersect_file from $RefSeq_CDS_FILE and $TARGET_FILE"
/hgsc_software/BEDTools/latest/bin/bedtools intersect -a $RefSeq_CDS_FILE -b $TARGET_FILE > $exon_target_intersect_file

