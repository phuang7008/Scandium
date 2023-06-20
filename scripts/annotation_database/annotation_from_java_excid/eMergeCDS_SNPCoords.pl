#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the RefSeq annotation\n";
my $type = shift || die "Please enter the version of annotation hg38 or hg37\n";

my $db_gene = $type=~/hg38/i ? "eMerge_CDS38" : "eMerge_CDS37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS $db_gene");
};

# now create table again
$dbh->do(qq{
  CREATE TABLE $db_gene (
  `id`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `cds_target_start` INT UNSIGNED NOT NULL,
  `cds_target_end` INT UNSIGNED NOT NULL,
  `exon_id` INT SIGNED NOT NULL,
  `exon_count` INT UNSIGNED NOT NULL,
  `cds_start` INT UNSIGNED NOT NULL,
  `cds_end` INT UNSIGNED NOT NULL,
  `cds_length` INT SIGNED NOT NULL,
  `gene_symbol`    varchar(200) NOT NULL,
  `gene_name`    varchar(200) NOT NULL,
  PRIMARY KEY (id),
  INDEX `CHR` (`chrom`),
  INDEX START (cds_start),
  INDEX END (cds_end),
  INDEX `COMP` (`chrom`, cds_target_start, cds_target_end),
  INDEX `SYM` (`gene_symbol`),
  INDEX `NAME` (`gene_name`)
) ENGINE=InnoDB;
});

my $prev_gene_name=".";
my $cds_length=0;
my $exon_count=0;

open (IN, $file) or die "open failed for reading: $!";

while (<IN>) {
	chomp;
	my @items = split "\t";

	# process the 3rd item on the list that contains the cds_id, gene name and gene symbol
	#
	my ($gene_symbol, $cds_info) = split(/\|/, $items[3]);
	my ($gene_name, $cds_id);

	if ($cds_info=~/cds/i) {
		($gene_name, $cds_id) = split(/_cds_/, $cds_info);
	} else {
		($gene_name, $cds_id) = ($cds_info, 0);
	}

	if ($gene_name eq $prev_gene_name) {
		# refseq/CCDS/vega gene is 0 based (bed file format), so don't need to + 1
		# but it seems the old ExCID + 1, so I will add 1 here
		#
		$cds_length += $items[2] - $items[1] + 1;
		$exon_count += 1;
	} else {
		if ($prev_gene_name ne ".") {
			if ($prev_gene_name=~/^NM/) {
				# only need to update the gene's cds_length info for the corresponding rows
				#
				$sql = "UPDATE $db_gene SET exon_count=$exon_count, cds_length=$cds_length  WHERE gene_name='$prev_gene_name'";
				$sth = $dbh->prepare($sql) or die "Query problem $!\n";
				$sth->execute() or die "Execution problem $!\n";
			}
		}

		# start a new gene cds
		#
		$cds_length = $items[2] - $items[1] + 1;
		$exon_count = 1;
		$prev_gene_name = $gene_name;
	}

	# now do the insertion
	#
	if ($gene_name=~/^NM/) {
		$sql = "INSERT INTO $db_gene VALUES (0, '$items[0]', $items[1], $items[2], $cds_id, $exon_count, 0, 0, 0, '$gene_symbol', '$gene_name')";
		$sth = $dbh->prepare($sql) or die "Query problem $!\n";
		$sth->execute() or die "Execution problem $!\n";
	} else {
		# for snp
		#
		$sql = "INSERT INTO $db_gene VALUES (0, '$items[0]', $items[1], $items[2], -1, 0, $items[1]-1, $items[2], 1, '$gene_symbol', '$gene_name')";
		$sth = $dbh->prepare($sql) or die "Query problem $!\n";
		$sth->execute() or die "Execution problem $!\n";
	}
}
close(IN);

# the last entry is a SNP, so we are fine!

print("Finish DB dumping\n");

$dbh->disconnect();

exit;
