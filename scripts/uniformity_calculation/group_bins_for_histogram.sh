#!/bin/bash

# this script is used to group bin count together for histogram plotting!
#
if [[ $# -ne 0 ]]; then
	echo "Illegal Number of Parameters"
	echo "/stornext/snfs5/next-gen/scratch/phuang/projects/group_related_works/Novaseq_bias/group_bins_for_histogram.sh"
	exit
fi

path=`pwd`
resultPath=$path/binned_groups
echo "$path"

id=1
lower_bound=2
upper_bound=150

for bFile in `ls $path/*WGS_between*_REPORT.txt`
#for bFile in `ls $path/*wo_centromeres`
do
	bambasename=""
	if [[ $bFile = *"_wo_centromeres"* ]]; then                                                               
		bambasename=$(basename "$bFile" _wo_centromeres)                              
	fi

	if [[ $bFile = *"hgv.cram"* ]]; then
		bambasename=$(basename "$bFile" .hgv.cram.WGS_between1x_150x_REPORT.txt) 
	fi

	if [[ $bFile = *"hgv.bam"* ]]; then 
		bambasename=$(basename "$bFile" .hgv.bam.WGS_between1x_150x_REPORT.txt) 
	fi

	if [[ $bFile = *"realigned.recal.bam"* ]]; then
		bambasename=$(basename "$bFile" .realigned.recal.bam.WGS_between1x_150x_REPORT.txt) 
	fi

	if [[ $bFile = *"realigned.recal.cram"* ]]; then
		bambasename=$(basename "$bFile" .realigned.recal.cram.WGS_between1x_150x_REPORT.txt) 
	fi

	#echo "$bambasename"

	outfile=${bambasename}_binned
	#echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/group_bins_together.pl $bFile > $resultPath/${bambasename}_binned"
	echo "/hgsc_software/perl/perl-5.18.2/bin/perl /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/uniformity_calculation/group_bins_together.pl $bFile > $resultPath/$outfile" | msub -V -d $resultPath -q high -A proj-dm0001 -j oe -N group_bins-$id -l nodes=1:ppn=1,mem=2gb

	# now we need to draw the histogram here
	echo "/hgsc_software/R/R-3.4.3/bin/Rscript /stornext/snfs5/next-gen/scratch/phuang/git_repo/graph_scripts/R_scripts/distribution_coverage_binning_smoothed.R $resultPath/$outfile $lower_bound $upper_bound" | msub -V -d $resultPath -q high -A proj-dm0001 -j oe -N histogram-$id -l depend=afterok:group_bins-$id -l nodes=1:ppn=1,mem=8gb

	id=$((id+1))

	#break
done
