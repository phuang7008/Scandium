#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the file name contains the RefSeq CDS annotation\n";
my $type = shift || die "Please enter the version of annotation hg38 or hg37\n";
my $database = $type=~m/hg38/i ? "Gene_RefSeq_CDS38" : "Gene_RefSeq_CDS37";
my $hgnc = $type eq "hg38" ? "HGNC38" : "HGNC37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS $database");
};

# now create table again
$dbh->do(qq{
  CREATE TABLE $database (
  `id`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `exon_start` INT UNSIGNED NOT NULL,
  `exon_end`   INT UNSIGNED NOT NULL,
  `annotation` varchar(4200) NOT NULL,
  PRIMARY KEY (id),
  INDEX `CHR` (`chrom`),
  INDEX `COMP` (`chrom`, exon_start, exon_end),
  INDEX `ANN` (`annotation`)
) ENGINE=InnoDB;
});

my ($sql, $sth);

open (IN, $file) or die "open failed for reading: $!";
#my $header=<IN>;

while (<IN>) {
	chomp;
	my @items = split "\t";

	# now do the insertion
	$sql = "INSERT INTO $database VALUES (0, '$items[0]', $items[1], $items[2], '$items[3]')";
	$sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";
}

print("Finish DB dumping\n");

$dbh->disconnect();
close(IN);

exit;
