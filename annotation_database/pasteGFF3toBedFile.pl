#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;

my $file = shift || die "Please enter the name of the file that contains all the intronic regions\n";

open (IN, $file) or die "open failed for reading: $!";

while (<IN>) {
	next if !/^chr/;
	chomp;
	my @items = split "\t";

	# 0=chr; 3=start; 4=end; 8=gene desc
	my @desc = split(/;/, $items[8]);
	foreach my $dc (@desc) {
		if ($dc=~/Name=(.*)/) {
			my $name = $1;
			my $alias = uc($name);
			$alias=~s/HSA-MIR-/MIR/;
			#$alias=~tr/MIR_/MIR/;
			print("$items[0]\t$items[3]\t$items[4]\t$name\t$alias\n");
		}
	}
}

close(IN);
