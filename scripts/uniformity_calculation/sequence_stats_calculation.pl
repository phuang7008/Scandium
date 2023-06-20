#!/hgsc_software/perl/perl-5.18.2/bin/perl
# This script is used to calculate the total count for each bin
# from between_lowX_highX file
#
use strict;
use Data::Dumper;
use File::Basename;
use POSIX;

my $file = shift || die "please specify the input file name\n";
#my $ave_in = shift || die "please specify the average coverage for the run\n";
my $out_lier_cutoff = shift;
$out_lier_cutoff = 1000 if !defined $out_lier_cutoff;

my $bname = basename($file);
my $ave = &average($file);
my ($std, $var) = &stdev($file, $ave);

# calculate coefficient of variation
#
my $cv = $std / $ave;

#print("sample\tMean\tstdev\tabs_variance\tcoefficient_of_variation\n");
print("$bname\t$ave\t$std\t$var\t$cv\n");
#print("The raw data average is $ave_in\n");
#print("The smoothed data average is $ave\n");
#print("The smoothed data std: $std and absolute variance is $var\n");

###################################################
sub average{
	my($file) = @_;
	open(IN, $file) || die "open failed for reading:$!";

	my $total = 0;
	my $count = 0;
	while (<IN>) {
		next if /^#/;
		chomp;
		my @vals = split "\t";

		# need to handle the out_liers
		#
		if ($vals[4] <= $out_lier_cutoff) {
			$total += $vals[3] * $vals[4];
			$count += $vals[3];
		}
	}
	close(IN);
        
	my $average = $total / $count;
	return $average;
}

sub stdev{
	my($file, $average) = @_;
	#my $average = &average($file);
	open(IN, $file) || die "open failed for reading:$!";

	my $sqr_total = 0;
	my $abs_total = 0;
	my $count = 0;
	while (<IN>) {
		next if /^#/;
		chomp;
		my @vals = split "\t";

		# need to handle out_liers
		#
		if ($vals[4] <= $out_lier_cutoff) {
			foreach my $v (1..$vals[3]) {
				$sqr_total += ($average-$vals[4]) ** 2;
				$abs_total += abs($average-$vals[4]);
			}
			$count += $vals[3];
		}
	}

	my $std = ($sqr_total / $count) ** 0.5;
	my $var = ($abs_total / $count);

	return ($std, $var);
}
