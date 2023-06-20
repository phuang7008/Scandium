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

$version = 'hg38' if !defined $version;
$a_size  = 7 if !defined $a_size;

# need to remove Ns regions for the calculation
#
my $num_of_Ns;

if ($version eq "hg38") {
	$num_of_Ns = 173893331;
} else {
	$num_of_Ns = 237019493;
}

my %hist;
my $total_bases = 0;
open (IN, $file_in) or die "open $file_in failed$!";

while (<IN>) {
	next if /^#/;

	chomp;
	my ($chr_id, $start, $stop, $len, $cov) = split "\t";
	$hist{$cov} += $len;
	$total_bases += $len;
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
$baseFile=~s/\.hgv\.cram\.WGS\_between1x\_150x\_REPORT\.txt//;
$baseFile=~s/\.hgv\.bam\.WGS\_between1x\_150x\_REPORT\.txt//;
$baseFile=~s/\.realigned\.recal\.bam\.WGS\_between1x\_150x\_REPORT\.txt//;
$baseFile=~s/\.realigned\.recal\.cram\.WGS\_between1x\_150x\_REPORT\.txt//;

# now need to remove the Ns regions
#
#my $pct = $mode_auc / $total_auc;		# include Ns
#$total_auc -= $num_of_Ns;

my $p_pct = int (1000 * ($mode_auc / $total_auc) + 5);
$p_pct = $p_pct / 1000;

my $m_pct = int (1000 * ($median_auc / $total_auc) + 5);
$m_pct = $m_pct / 1000;

my $mean_pct = int (1000 * ($mean_auc / $total_auc) + 5);
$mean_pct = $mean_pct / 1000;

# now we need to calculate kurtosis
#
my $mean = 0;
my $std = 0;
my $n = 0;
foreach my $cov (keys %hist) {
	foreach my $d (1..$hist{$cov}) {
		$mean += $cov;
		$std += $cov*$cov;
		$n += 1;
	}
}
$mean = $mean / $n;
$std = $std / $n - $mean*$mean; 
$std = sqrt($std); 

my $x = 0;
foreach my $cov (keys %hist) {
	foreach my $d (1..$hist{$cov}) {
		$x += pow(($cov - $mean) / $std, 4);		
	}
}
my $kurt = $x/$n - 3;

# using different approach
#
my $m2 = 0;
my $m3 = 0;
my $m4 = 0;

foreach my $cov (keys %hist) {
    foreach my $d (1..$hist{$cov}) {
        $m2 += pow($cov - $mean, 2);
        $m3 += pow($cov - $mean, 3);
        $m4 += pow($cov - $mean, 4);
    }
}

$m2 /= $total_bases;   
$m3 /= $total_bases;
$m4 /= $total_bases;

# population kurtosis
#
my $kurt_p = ($m4/ pow($m2, 2)) -3;  

# for sample kurtosis
#
my $kurt_s = (($total_bases+1)*$kurt_p)+6;   
$kurt_s *= ($total_bases-1)/(($total_bases-2)*($total_bases-3));

# output input parameters for reading to review
#
print("filename\tref_version\tsize_picked\tMean\tUniformity_Mean\tMode(Peak)\tUniformity_Mode\tMedian\tUniformity_Median\tKurtosis(Pop)\tKurtosis(Sample)\tKurtosis(Z)\n");
print("$baseFile\t$version\t$a_size\t$mean\t$mean_pct\t$mode_cov_idx\t$p_pct\t$median\t$m_pct\t$kurt_p\t$kurt_s\t$kurt\n");

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
