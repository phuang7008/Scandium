#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;

# in order to draw the graph on R, I need to create a matrix.
# in order to create a matrix, I need to know the dimention
# here I will try to find out the length needed to draw the graph

my $OMIM_transcript_sorted_bed = shift || die "Please specify the OMIM transcript sorted bed file here\n";
my $chrom_id = shift || die "please enter the chromosome id you are working on\n";

open (IN, $OMIM_transcript_sorted_bed) or die "open failed for reading: $!";

my $total_length;

while (<IN>) {
	chomp;
	my @items = split "\t";
	next if $items[0] ne $chrom_id;
	#next if length(@items) < 6;

	$total_length = $total_length + $items[2] - $items[1] + 1;
	#$total_length = $total_length + $items[3] + 1;
	#$total_length = $total_length + $items[3];
}

print $total_length;

close(IN);

exit;
