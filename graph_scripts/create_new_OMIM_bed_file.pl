#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;

# this script is used to generate a bed file for all OMIM genes in transcript format (not exon format)
# Therefore, the script will just eliminate all the exon info

my $in_file  = shift || die "Please specify the OMIM database file\n";
my $out_file = shift || die "Please specify the output file name for OMIM transcript bed file\n";

#open (IN, "/hgsc_software/production/users/cbuhay/ExCID/latest/database/PKv2_Gene_Database.bed") or die "open failed for reading: $!";
open (IN, $in_file) or die "open failed for reading: $!";

# for %genes hash, there are 4 levels: 
# first key: chrom id
# second key: gene symbol
# third key: transcript name
# fourth key: transcript start position
# fifth key: transcript end position
#
my %genes;

while (<IN>) {
	chomp;
	my @items = split "\t";
	my ($gene_symbol, $gene_name) = split '\|', $items[3];
	$gene_name=~s/_cds_\d+//;

	if (exists $genes{$items[0]} && $genes{$items[0]}{$gene_symbol} && $genes{$items[0]}{$gene_symbol}{$gene_name}) {
		if ($genes{$items[0]}{$gene_symbol}{$gene_name}{start} > $items[1]) {
			$genes{$items[0]}{$gene_symbol}{$gene_name}{start} = $items[1];
		}

		if ($genes{$items[0]}{$gene_symbol}{$gene_name}{end} < $items[2]) {
			$genes{$items[0]}{$gene_symbol}{$gene_name}{end} = $items[2];
		}
	} else {
		$genes{$items[0]}{$gene_symbol}{$gene_name}{start} = $items[1];
		$genes{$items[0]}{$gene_symbol}{$gene_name}{end} = $items[2];
	}
}

close(IN);

# now output the hashtable to a bed file
#
open (OUT, ">$out_file") or die "open failed for writing: $!";

foreach my $cid (sort keys %genes) {
	foreach my $sym (keys %{$genes{$cid}}) {
		foreach my $name (keys %{$genes{$cid}{$sym}}) {
			print OUT "$cid\t$genes{$cid}{$sym}{$name}{start}\t$genes{$cid}{$sym}{$name}{end}\t$sym-$name\n";
		}
	}
}

close(OUT);

exit;
