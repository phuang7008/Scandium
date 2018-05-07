#!/hgsc_software/perl/perl-5.18.2/bin/perl
# This script is used to calculate the total count for each bin
# from between_lowX_highX file
# It finds the max value as the center of the peak and get the Peak area under the curve
# based on user's size input. The peak size is the total points used for the calculation
# Therefore, in this case, the peak doesn't have to be symmetric!
#
use strict;
use Data::Dumper;
use File::Basename;
use POSIX;

my $file = shift || die "please specify the input file name\n";
my $peak_size = shift || die "please specify the number of points around the Peak that would be used for the Peak AUC calculation\n";
my $version   = shift || die "Please specify the version of reference used (hg19 or hg38)\n";
my $num_of_Ns;

if ($version eq "hg38") {
	$num_of_Ns = 173893331;
} else {
	$num_of_Ns = 237019493;
}

my %hist;
open (IN, $file) or die "open $file failed$!";

while (<IN>) {
	next if /^#/;

	chomp;
	my ($chr_id, $start, $stop, $len, $cov) = split "\t";
	$hist{$cov} += $len;
}
close IN;

# find the peak, that is the max
my $peak_cov_count=0;
my $peak_cov_idx=0;
foreach my $cov (sort {$a <=> $b} keys %hist) {
	#print("$cov\t$hist{$cov}\n");
	next if $cov == 0;

	if ($peak_cov_count < $hist{$cov}) {
		$peak_cov_count = $hist{$cov};
		$peak_cov_idx   = $cov;
	}
}

# calculate area under curve
#
my ($total_auc, $peak_auc) = (0,0);

foreach my $cov (sort {$a <=> $b} keys %hist) {
	$total_auc += $hist{$cov};	
}

$peak_auc = $peak_cov_count;
my $counter=1;
my $left  = $peak_cov_idx-1;
my $right = $peak_cov_idx+1;
my @peak_selected=();

while (1) {
	if ($hist{$left} > $hist{$right}) {
		$peak_auc += $hist{$left};
		push @peak_selected, $hist{$left};

		$counter += 1;
		$left -= 1;
	} elsif ($hist{$left} < $hist{$right}) {
		$peak_auc += $hist{$right};
		push @peak_selected, $hist{$right};

		$counter += 1;
		$right += 1;
	} else {
		if ($counter - $peak_size >= 2) {
			$peak_auc += $hist{$left} * 2;
			push @peak_selected, $hist{$left};
			push @peak_selected, $hist{$left};

			$left -= 1;
			$right += 1;
			$counter += 2;
		} else {
			$peak_auc += $hist{$left};
			push @peak_selected, $hist{$left};

			$counter += 1;
			$left -= 1;
		}
	}

	if ($counter == $peak_size) {
		last;
	}
}

# get the file name
#
my $baseFile = basename($file);
$baseFile=~s/\.hgv\.cram\.WGS\_between1x\_150x\_REPORT\.txt//;
$baseFile=~s/\.hgv\.bam\.WGS\_between1x\_150x\_REPORT\.txt//;
$baseFile=~s/\.realigned\.recal\.bam\.WGS\_between1x\_150x\_REPORT\.txt//;
$baseFile=~s/\.realigned\.recal\.cram\.WGS\_between1x\_150x\_REPORT\.txt//;

# now need to remove the Ns regions
#
#print("Total area under the curve with Ns is $total_auc\n");
my $pct = $peak_auc / $total_auc;
#print("The percentage include Ns is $pct\n");

$total_auc -= $num_of_Ns;

$pct = int (100000 * ($peak_auc / $total_auc));
$pct = $pct / 1000;

print("$baseFile\t$pct\n");
#print("total auc without Ns is $total_auc\n");
#print("peak auc is $peak_auc\n");
#print("The percentage without Ns is $pct\n");

#print("The peak point is $peak_cov_idx\n");
#print("The peak coverage count is $peak_cov_count\n");
#print("The peak range selected:\n");

#foreach my $val (@peak_selected) {
#	print("\t$val\n");
#}
