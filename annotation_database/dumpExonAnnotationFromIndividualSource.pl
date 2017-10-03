#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
# this script is used to dump the exon annotations from different sources 
# such as RefSeq, CCDS, VEGA or Gencode etc.
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file   = shift || die "Please enter the name of the file that contains merged exon annotation\n";
my $type   = shift || die "Please specify the version of gene annotation hg19 or hg38\n";
my $source = shift || die "Please specify the source of the annotations such as RefSeq, CCDS, VEGA or Gencode etc.\n";
if ($source=~/refseq/i) {
	$source = "RefSeq";
} elsif ($source=~/CCDS/i) {
	$source = "CCDS";
} elsif ($source=~/vega/i) {
	$source = "VEGA";
} elsif ($source=~/gencode/i) {
	$source = "Gencode";
} else {
	print "The source $source is not defined\n";
	exit;
}

my $database = $type eq "hg38" ? $source."_Exon38" : $source."_Exon37";
my $hgnc = $type eq "hg38" ? "HGNC38" : "HGNC37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS $database");
};

# now create table again
$dbh->do(qq{
CREATE TABLE $database (
  `eid`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `exon_start`  INT UNSIGNED NOT NULL,
  `exon_end`    INT UNSIGNED NOT NULL,
  `annotation` varchar(3000) NOT NULL,
  PRIMARY KEY (eid),
  INDEX `COMP` (`chrom`, exon_start, exon_end),
  INDEX `ANNO` (`annotation`)
) ENGINE=InnoDB;
});

open (IN, $file) or die "open failed for reading: $!";

while (<IN>) {
	chomp;
	my @items = split "\t";

	my $sql = "INSERT INTO $database VALUES (0, '$items[0]', $items[1], $items[2], '$items[3]')";
	my $sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";
}

$dbh->disconnect();
close(IN);

exit;
