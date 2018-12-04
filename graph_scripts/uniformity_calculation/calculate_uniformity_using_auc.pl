#!/hgsc_software/perl/perl-5.18.2/bin/perl
# This script is used to calculate the total count for each bin
# from between_lowX_highX file
#
use strict;
use Data::Dumper;
use POSIX;

my $file = shift || die "please specify the input uniformity file name\n";
my $average = shift || die "please specify the average coverage for the run\n";
my $version = shift || die "Please specify the version of reference used (hg19 or hg38)\n";
my $num_of_Ns;

if ($version eq "hg38") {
	$num_of_Ns = 173893331;
} else {
	$num_of_Ns = 237019493;
}

#my $upper_bound = ceil ($average + 3);
my $upper_bound = ceil ($average + 4);
#my $upper_bound = ceil ($average + 5);
#my $lower_bound = floor ($average - 5);
my $lower_bound = floor ($average - 4);
#my $lower_bound = floor ($average - 3);

my %hist;
open (IN, $file) or die "open $file failed$!";

while (<IN>) {
	next if /^#/;

	chomp;
	my ($chr_id, $start, $stop, $len, $cov) = split "\t";
	$hist{$cov} += $len;
}
close IN;

foreach my $cov (sort {$a <=> $b} keys %hist) {
	print("$cov\t$hist{$cov}\n");
}

# calculate area under curve
#
my ($total_auc, $peak_auc) = (0,0);

foreach my $cov (sort {$a <=> $b} keys %hist) {
	#$total_auc += 0.5 * ($cov - $prev_cov) * ($prev_len + $hist{$cov});	
	$total_auc += $hist{$cov};	

	if ( $lower_bound <= $cov && $cov <= $upper_bound) {
		#$peak_auc += 0.5 * 1 * ($prev_len + $hist{$cov});
		$peak_auc += $hist{$cov};
	}
}

# now need to remove the Ns regions
#
print("Total area under the curve with Ns is $total_auc\n");
my $pct = $peak_auc / $total_auc;
print("The percentage include Ns is $pct\n");

$total_auc -= $num_of_Ns;

$pct = $peak_auc / $total_auc;
print("total auc without Ns is $total_auc\n");
print("peak auc is $peak_auc\n");
print("The percentage without Ns is $pct\n");
print("lower bound is $lower_bound\n");
print("upper bound is $upper_bound\n");
