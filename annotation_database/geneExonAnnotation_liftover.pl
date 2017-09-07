#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the RefSeq annotation\n";
my $type = shift || die "Please enter the type of the analysis hg19 vs hg38\n";
my $database = $type eq "hg38" ? "Capture_Regions38" : "Capture_Regions37";
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
  `id`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `start` INT UNSIGNED NOT NULL,
  `end` INT UNSIGNED NOT NULL,
  `exon_start` INT UNSIGNED NOT NULL,
  `exon_end` INT UNSIGNED NOT NULL,
  `target_start` INT UNSIGNED NOT NULL,
  `target_end` INT UNSIGNED NOT NULL,
  `exon_id` INT UNSIGNED NOT NULL,
  `exon_count` INT UNSIGNED NOT NULL,
  `gene_symbol`    varchar(200) NOT NULL,
  `refseq_name`    varchar(200) NOT NULL,
  PRIMARY KEY (id),
  INDEX `CHR` (`chrom`),
  INDEX START (target_start),
  INDEX END (target_end),
  INDEX EXON_START (exon_start),
  INDEX EXON_END (exon_end),
  INDEX `COMP` (`chrom`, start, end),
  INDEX `SYM` (`gene_symbol`),
  INDEX `NAME` (`refseq_name`)
) ENGINE=InnoDB;
});

open (IN, $file) or die "open failed for reading: $!";
#my $header=<IN>;

while (<IN>) {
	chomp;
	my @items = split "\t";

	# process the 3rd item on the list that contains the exon_id, refseq name and gene symbol
	my @info = split(/;/, $items[3]);

	foreach my $region (@info) {
		my ($ref_exon_id, $gene_symbol, $exon_start, $exon_end) = split(/=/, $region);
		my ($refseq_name, $exon_id);
	   	if ($ref_exon_id=~/(^N.*)_(\d+)$/) {
			($refseq_name, $exon_id) = ($1, $2);
		} elsif ($ref_exon_id=~/(^X.*)_(\d+)$/) {
			($refseq_name, $exon_id) = ($1, $2);
		}
		
		my ($sql, $sth);
		if (!defined $gene_symbol) {
			$sql = "SELECT symbol FROM $hgnc WHERE (refseq_accession IS NOT NULL AND find_in_set('$refseq_name', refseq_accession))";
			$sth = $dbh->prepare($sql) or die "DB query error: $!";
			$sth->execute() or die "DB execution error: $!";

			($gene_symbol) = $sth->fetchrow_array if (!defined $gene_symbol);
		}

		# now do the insertion
		$sql = "INSERT INTO $database VALUES (0, '$items[0]', $items[1], $items[2], $exon_start, $exon_end, $items[5], $items[6], $exon_id, 1, '$gene_symbol', '$refseq_name')";
		$sth = $dbh->prepare($sql) or die "Query problem $!\n";
		$sth->execute() or die "Execution problem $!\n";
	}
}

print("Finish DB dumping. Now update exon count info\n");

# now I need to update the exon_count information
$sql = "select refseq_name, count(exon_id) from $database group by refseq_name";
$sth=$dbh->prepare($sql) or die "DB query error: $!";
$sth->execute() or die "DB execution error: $!";

my %refseq;
while ( my($refseq_name, $count) = $sth->fetchrow_array) {
    $refseq{$refseq_name} = $count;
}

foreach my $refseq_name (keys %refseq) {
    # now update the table
    $sql = "UPDATE $database SET exon_count=$refseq{$refseq_name} WHERE refseq_name='$refseq_name'";
    $sth = $dbh->prepare($sql) or die "Query problem $!\n";
    $sth->execute() or die "Execution problem $!\n";
}

$dbh->disconnect();
close(IN);

exit;
