#!/usr/bin/perl
#
# some refseq sequence names are the same even though they point to different genomic location
# to eliminate the problem, we need to add the suffix to separate them apart
# we will eventually remove these suffix in another script (need to be handled later!)
#
use FindBin;
use lib "$FindBin::Bin";
use strict;

my $file = shift || die "Please enter the name of the refseq file\n";

open (IN, $file) || die "open file $file failed: $!";
my $header=<IN>;

my %refseq;
while (<IN>) {
	my @items = split /\t/;
	$refseq{$items[0]}++
}
close(IN);

open(IN, $file) || die "open file $file failed: $!";
$header=<IN>;

my %new_refseq;
while (<IN>) {
	my @items = split /\t/;

	if ($refseq{$items[0]} >= 2) {
		my $suffix = $new_refseq{$items[0]} + 1;
		$new_refseq{$items[0]}++;
		$items[0] .= "-$suffix";
	} else {
		$new_refseq{$items[0]}++;
	}

	my $line = join("\t", @items);
	print $line;
}
close(IN);

exit;
