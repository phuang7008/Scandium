#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the MANE file $!\n";
my $database = "MANE_GRCh38";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

# drop the table first
eval {
    $dbh->do("DROP TABLE IF EXISTS $database");
};

# now create table again
$dbh->do(qq{
CREATE TABLE $database (
  `mane_id` INT UNSIGNED NOT NULL,
  `symbol`  varchar(150) NOT NULL,
  `alias_symbol` varchar(150) NULL,
  `refseq_transcript` varchar(150) NULL,
  `ensembl_transcript` varchar(150) NULL,
  PRIMARY KEY (mane_id),
  INDEX `REF` (`refseq_transcript`),
  INDEX `EMSEMBL` (`ensembl_transcript`),
  INDEX `SYM` (`symbol`)
) ENGINE=MyISAM;
});

open(IN, $file) or die "read failed for $file: $!";
my $sql;
my $sth;
my %seen;
my $header = <IN>;

while(<IN>) {
	chomp;
	my @items = split("\t");
    next if $items[0] eq "";
    next if $seen{$items[0]};

	if ($items[1] eq "") { $items[1]=undef; }
	if ($items[2] eq "") { $items[2]=undef; }
	if ($items[3] eq "") { $items[3]=undef; }
	if ($items[4] eq "") { $items[4]=undef; }
    $seen{$items[0]}++;

	$sql = "INSERT INTO $database (mane_id, symbol, alias_symbol, refseq_transcript, ensembl_transcript) VALUES (?, ?, ?, ?, ?)";
	$sth = $dbh->prepare($sql) or die "DB query error: $!";
	$sth->execute($items[0], $items[2], $items[1], $items[3], $items[4]) or die "execution failed: $dbh->errstr()";
}

$dbh->disconnect();
close(IN);

exit;
