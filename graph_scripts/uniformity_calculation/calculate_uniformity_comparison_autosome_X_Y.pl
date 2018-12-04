#!/hgsc_software/perl/perl-5.18.2/bin/perl
# this script is used to calculation the uniformity metric using both median and mode (mode)
#
# the median value will have to be part of the input value
#
# For the mode of the mode:
# It finds the max value as the center of the mode and get the Peak area under the curve
# based on user's size input. The mode size is the total points used for the calculation
# Therefore, in this case, the mode doesn't have to be symmetric!
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/';

use strict;
use Data::Dumper;
use File::Basename;
use POSIX;
use Statistics::Basic qw(:all);

my $file_in = shift || die "please specify the input uniformity file name\n";
my $version = shift;
my $a_size  = shift;	# the mode area size
my $XY = shift;			# include sex chromosome or not (1 for yes, 2 for no)

$version = 'hg38' if !defined $version;
$a_size  = 7 if !defined $a_size;

# need to remove Ns regions for the calculation
#
my $num_of_Ns = 237019493;	# for hg37

if ($version eq "hg38") {
	$num_of_Ns = 173893331;
}

if ($XY != 1) {
	# need to modify Ns regions (not include X and Y)
	#
	if ($version eq "hg38") {
		$num_of_Ns -= 34738921;
	} else {
		$num_of_Ns -= 40559035;
	}
}

my @chrom_list=("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
	              "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT");

if ($version eq "hg38") {
	@chrom_list=("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
		"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
		"chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM");
}

my %chrom_hash;
foreach my $id (@chrom_list) {
	$chrom_hash{$id}++;
}

my %hist;
my $total_bases = 0;
open (IN, $file_in) or die "open $file_in failed$!";

while (<IN>) {
	next if /^#/;

	chomp;
	my ($chr_id, $start, $stop, $len, $cov) = split "\t";
	next if $XY != 1 && ($chr_id eq "X" || $chr_id eq "Y" || $chr_id eq "chrX" || $chr_id eq "chrY") ;

	if (defined $chrom_hash{$chr_id}) {
		$hist{$cov} += $len;
		$total_bases += $len;
	}
}
close IN;

# need to remove Ns regions
#
$total_bases = $total_bases - $num_of_Ns;
$hist{0} = $hist{0} - $num_of_Ns;

my $half_total = 0;
my $median = 0;
foreach my $key (sort {$a <=> $b} keys %hist) {
	if ($half_total == ($total_bases/2)) {
		$median = $key;
		last;
	} elsif ($half_total > ($total_bases/2)) {
		$median = $key - 1;
		last;
	} else {
		$half_total += $hist{$key};
	}
}

# find the mode (Mode), that is the max
#
my $mode_cov_count=0;
my $mode_cov_idx=0;
my $sum=0;
foreach my $cov (sort {$a <=> $b} keys %hist) {
	if ($cov == 0) {
		#$sum += $hist{$cov};
		next;
	}

	if ($mode_cov_count < $hist{$cov}) {
		$mode_cov_count = $hist{$cov};
		$mode_cov_idx   = $cov;
	}

	$sum += $hist{$cov} * $cov;;
}

my $mean = int $sum / $total_bases + 0.5;

# calculate area under curve
#
my $total_auc = 0;

foreach my $cov (sort {$a <=> $b} keys %hist) {
	$total_auc += $hist{$cov};	
}

my $mean_auc   = &dynamicCalculateArea($mean);
my $mode_auc   = &dynamicCalculateArea($mode_cov_idx);
my $median_auc = &dynamicCalculateArea($median);

# get the file name
#
my $baseFile = basename($file_in);
$baseFile=~s/\.hgv\.cram\.WGS\_uniformity\_REPORT\.txt//;
$baseFile=~s/\.hgv\.bam\.WGS\_uniformity\_REPORT\.txt//;
$baseFile=~s/\.realigned\.recal\.bam\.WGS\_uniformity\_REPORT\.txt//;
$baseFile=~s/\.realigned\.recal\.cram\.WGS\_uniformity\_REPORT\.txt//;
$baseFile=~s/\.marked\.bam\.WGS\_uniformity\_REPORT\.txt//;

my $p_pct = int (1000 * ($mode_auc / $total_auc) + 0.5);
$p_pct = $p_pct / 1000;

my $m_pct = int (1000 * ($median_auc / $total_auc) + 0.5);
$m_pct = $m_pct / 1000;

my $mean_pct = int (1000 * ($mean_auc / $total_auc) + 0.5);
$mean_pct = $mean_pct / 1000;

# output input parameters for reading to review
#
print("filename\tref_version\tsize_picked\tMean\tUniformity_Mean\tMode(Peak)\tUniformity_Mode\tMedian\tUniformity_Median\n");
print("$baseFile\t$version\t$a_size\t$mean\t$mean_pct\t$mode_cov_idx\t$p_pct\t$median\t$m_pct\n");

# for debugging purpose
#
#print("Mode area under histogram is $mode_auc\n");
#print("Total area under histogram is  $total_auc\n");

#foreach my $cov (sort {$a <=> $b} keys %hist) {
#	print("$cov\t$hist{$cov}\n");
#}

exit;

####################################################################################
# Sub-routines
####################################################################################
sub dynamicCalculateArea {
	my ($cov_in) = @_;
	my $left  = $cov_in - 1;
	my $right = $cov_in + 1;

	my $cur_auc = $hist{$cov_in};
	my $counter = 1;

	if ($a_size == 1) {
		return $cur_auc;
	}

	while (1) {
		if ($hist{$left} > $hist{$right}) {
			$cur_auc += $hist{$left};
			$counter += 1;                                                                                    
			$left -= 1;
		} elsif ($hist{$left} < $hist{$right}) {
			$cur_auc += $hist{$right};
			$counter += 1;                                                                                    
			$right += 1;
		} else {
			# they are equal
			#
			if ($a_size - $counter >= 2) {
				$cur_auc += $hist{$left} * 2;
				$left -= 1;
				$right += 1;
				$counter += 2;
			} else {
				$cur_auc += $hist{$left};
				$left -= 1;
				$counter += 1;
			}
		}

		last if ($counter == $a_size);
	}

	return $cur_auc;
}
