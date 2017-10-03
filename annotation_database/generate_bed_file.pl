#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the file that contains all the intronic regions\n";
my $type = shift || die "Please specify the version of gene annotation hg19 or hg38\n";
my $hgnc = $type eq "hg38" ? "HGNC38" : "HGNC37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

open (IN, $file) or die "open failed for reading: $!";
while (<IN>) {
	next if /chrom/i;
	chomp;
	my @items = split "\t";

	my @starts = split ",", $items[6];
	my @ends   = split ",", $items[7];

	# remove 'chr' in front of chromosome id
	$items[1]=~s/chr//i;

	if ($items[0]=~/^CCDS.*/) {
		$items[0]=~s/(CCDS.*)\.\d+$/$1/;
		$sql = "SELECT symbol FROM $hgnc WHERE (ccds_id IS NOT NULL AND find_in_set('$items[0]', replace(ccds_id, '|', ',')))";
		$sth = $dbh->prepare($sql) or die "DB query error: $!";
		$sth->execute() or die "DB execution error: $!";

		my ($gene_symbol) = $sth->fetchrow_array;
		$items[8] = $gene_symbol;
		$sth->finish();
	}

	if (! defined $items[8]) { $items[8] = "."; }

	foreach my $idx (0..$#starts) {

		my $exon_id = $idx;
		if ($items[2] eq "-") {
			$exon_id = $items[5] - $idx - 1;
		}

		print "$items[1]\t$starts[$idx]\t$ends[$idx]\t$items[0]_exon_$exon_id\_$items[5]=$items[8]\t$items[8]\n";
	}
}

$dbh->disconnect();
close(IN);

exit;
