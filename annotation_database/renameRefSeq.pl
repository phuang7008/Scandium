#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
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
