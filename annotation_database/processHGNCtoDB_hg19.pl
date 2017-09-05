#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the HGNC file $!\n";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

# drop the table first
eval {
    $dbh->do("DROP TABLE IF EXISTS HGNC37");
};

# now create table again
$dbh->do(qq{
CREATE TABLE HGNC37 (
  `hgnc_id` INT UNSIGNED NOT NULL,
  `symbol`  varchar(150) NOT NULL,
  `alias_symbol` varchar(150) NULL,
  `prev_symbol` varchar(150) NULL,
  `vega_id` varchar(150) NULL,
  `ucsc_id` varchar(150) NULL,
  `refseq_accession` varchar(150) NULL,
  `ccds_id` varchar(150) NULL,
  `mirbase` varchar(150) NULL,
  PRIMARY KEY (hgnc_id),
  INDEX `REF` (`refseq_accession`),
  INDEX `CCDS` (`ccds_id`),
  INDEX `UCSC` (`ucsc_id`),
  INDEX `VEGA` (`vega_id`),
  INDEX `MIR` (`mirbase`),
  INDEX `SYM` (`symbol`)
) ENGINE=MyISAM;
});

open(IN, $file) or die "read failed for $file: $!";
my $sql;
my $sth;
my $header = <IN>;

while(<IN>) {
	chomp;
	my @items = split("\t");
	if ($items[1] eq "") { $items[1]=undef; }
	if ($items[2] eq "") { $items[2]=undef; }
	if ($items[3] eq "") { $items[3]=undef; }
	if ($items[4] eq "") { $items[4]=undef; }
	if ($items[5] eq "") { $items[5]=undef; }
	if ($items[6] eq "") { $items[6]=undef; }
	if ($items[7] eq "") { $items[7]=undef; }
	if ($items[8] eq "") { $items[8]=undef; }

	$sql = "INSERT INTO HGNC37 (hgnc_id, symbol, alias_symbol, prev_symbol, vega_id, ucsc_id, refseq_accession, ccds_id, mirbase) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)";
	$sth = $dbh->prepare($sql) or die "DB query error: $!";
	$sth->execute($items[0], $items[1], $items[2], $items[3], $items[4], $items[5], $items[6], $items[7], $items[8]);
}

$dbh->disconnect();
close(IN);

exit;
