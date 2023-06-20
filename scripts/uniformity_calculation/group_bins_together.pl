#!/hgsc_software/perl/perl-5.18.2/bin/perl
# This script is used to calculate the total count for each bin
# from between_lowX_highX file
#
my $file = shift || die "please specify the input file name\n";

my %hist;
open (IN, $file) or die "open $file failed$!";

while (<IN>) {
	next if /^#/;

	chomp;
	my ($chrom_id, $start, $end, $len, $cov) = split "\t";
	$hist{$cov} += $len;
}
close IN;

foreach $cov (sort {$a <=> $b} keys %hist) {
	print("$cov\t$hist{$cov}\n");
}

