#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the file with sorted and merged exons\n";
my $type = shift || die "Please enter the version of annotation hg38 or hg37\n";
my $database = $type=~/hg38/i ? "Gene_Exon38" : "Gene_Exon37";
my $hgnc = $type=~/hg38/i ? "HGNC38" : "HGNC37";
my $db2  = $type=~/hg38/i ? "Merged_Gene_Exon38" : "Merged_Gene_Exon37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS $database");
	$dbh->do("DROP TABLE IF EXISTS $db2");
};

# now create table again
$dbh->do(qq{
  CREATE TABLE $database (
  `ge_id`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `exon_start` INT UNSIGNED NOT NULL,
  `exon_end`   INT UNSIGNED NOT NULL,
  `exon_id`    INT UNSIGNED NOT NULL,
  `exon_count` INT UNSIGNED NOT NULL,
  `gene_symbol`  varchar(200) NOT NULL,
  `transcript_name`  varchar(200) NOT NULL,
  PRIMARY KEY (ge_id),
  INDEX `COMP` (`chrom`, exon_start, exon_end),
  INDEX `SYM`  (`gene_symbol`),
  INDEX `NAME` (`transcript_name`)
) ENGINE=InnoDB;
});

$dbh->do(qq{
  CREATE TABLE $db2 (
  `eid`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `exon_start` INT UNSIGNED NOT NULL,
  `exon_end`   INT UNSIGNED NOT NULL,
  `annotation`  varchar(3500) NOT NULL,
  PRIMARY KEY (eid),
  INDEX `COMP` (`chrom`, exon_start, exon_end),
  INDEX `ANN`  (`annotation`)
) ENGINE=InnoDB;
});

open (IN, $file) or die "open failed for reading: $!";
#my $header=<IN>;

while (<IN>) {
	chomp;
	my @items = split "\t";

	# first dump everything to db2
	my ($sql, $sth);
	$sql = "INSERT INTO $db2 VALUES (0, '$items[0]', $items[1], $items[2], '$items[3]')";
	$sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";
	
	undef $sql;
	undef $sth;

	# process the 3rd item on the list that contains the exon_id, refseq name and gene symbol
    # NR_046018|exon_0|gene_3=DDX11L1=11873=14409=1652
    #
	my @info = split(/;/, $items[3]);

	foreach my $region (@info) {
		my ($exon, $gene_symbol) = split(/=/, $region);
		my ($gene_name, $exon_id, $exon_count) = split(/\|/, $exon);
        $exon_id =~s/exon_//;
        if ($exon_count=~/gene/) {
            $exon_count =~s/gene_//;
        } elsif ($exon_count=~/miRNA/) {
            $exon_count =~s/miRNA_//;
        } else {
            print "Wrong annotation $exon_count\n";
        }

		
		# only CCDS doesn't have gene_symbol, we need to find it here!
		if (!defined $gene_symbol) {
			$gene_name=~s/(CCDS\d+)\.\d+/$1/;

			$sql = "SELECT symbol FROM $hgnc WHERE (ccds_id IS NOT NULL AND find_in_set('$gene_name', ccds_id))";
			$sth = $dbh->prepare($sql) or die "DB query error: $!";
			$sth->execute() or die "DB execution error: $!";

			($gene_symbol) = $sth->fetchrow_array if (!defined $gene_symbol);
		}

		# now do the insertion
		$sql = "INSERT INTO $database VALUES (0, '$items[0]', $items[1], $items[2], $exon_id, $exon_count, '$gene_symbol', '$gene_name')";
		$sth = $dbh->prepare($sql) or die "Query problem $!\n";
		$sth->execute() or die "Execution problem $!\n";
	}
}

print("Finish DB dumping.\n");

$dbh->disconnect();
close(IN);

exit;
