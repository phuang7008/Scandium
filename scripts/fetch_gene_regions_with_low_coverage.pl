#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;

my $low_cov_report_file = shift || die "Please specify the coverage report file\n";
#my $gene_list_file = shift || die "Please enter the file name contains OMIM gene list\n";
#open (IN, $gene_list_file) or die "open failed for reading: $!";

open (IN, "/hgsc_software/production/users/cbuhay/ExCID/latest/database/PKv2_Gene_Database.bed") or die "open failed for reading: $!";

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

# now find out the regions overlap
open (IN, $low_cov_report_file) or die "open failed for reading: $!";
my $header=<IN>;
$header=<IN>;
$header=<IN>;

while (<IN>) {
	chomp;
	my @items = split "\t";

	# process the 3rd item on the list that contains the exon_id, refseq name and gene symbol
	my @info = split(/;/, $items[3]);
}

close(IN);

exit;
