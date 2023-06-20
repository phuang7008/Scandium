#!/hgsc_software/perl/perl-5.18.2/bin/perl
# This script is used to calculate the total count for each bin
# from between_lowX_highX file
#
use strict;
use Data::Dumper;
use File::Basename;
use Statistics::Descriptive;

my $dir = shift || "Please specify the directory where the input files reside\n";
my $file_suffix = shift || die "please specify the input file suffix\n";

my @files = <"$dir/*$file_suffix">;

foreach my $file (@files) {
	my @my_data;
	open (IN, $file) or die "open $file failed$!";

	my $fname = basename($file);
    $fname=~s/\_distribution.txt$//;

	while (<IN>) {
		chomp;
		my ($cov, $count) = split "\t";
		foreach my $d (1..$count) {
			push @my_data, $cov;
		}
	}
	close IN;

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@my_data);

	open (OUT, ">$fname"."_stats.txt") or die "open failed for writing:$!";
	print OUT "Sample\tMean\tMode\tMedian\tKurtosis\n";
	print OUT "$fname\t";
	print OUT $stat->mean();
	print OUT "\t";
	print OUT $stat->mode();
	print OUT "\t";
	print OUT $stat->median();
	print OUT "\t";
	print OUT $stat->kurtosis();
	print "\n";

}

exit;
