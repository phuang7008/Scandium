#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the RefSeq annotation\n";
my $type = shift || die "Please enter the version of annotation hg38 or hg37\n";

my $db_gene = $type=~/hg38/i ? "PKv2_CDS38" : "PKv2_CDS37";

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
  `exon_id` INT UNSIGNED NOT NULL,
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

# because the file is sorted by the gene_symbol and gene_name
#
my $prev_gene_name=".";
my $cds_start=0;
my $cds_end=0;
my $cds_length=0;
my $exon_count=0;
my $gene_name="";

open (IN, $file) or die "open failed for reading: $!";

while (<IN>) {
	chomp;
	my @items = split "\t";

	# process the 4th item (index 3) on the list that contains the cds_id, refseq name and gene_symbol information
	# 22      43088895        43089957        A4GALT|NM_017436_cds_0
	#
	my ($gene_symbol, $refseq_info) = split(/\|/, $items[3]);
	my ($gene_name, $cds_id);
   	if ($refseq_info=~/cds/) {
	   	($gene_name, $cds_id) = split("_cds_", $refseq_info);
	} else {
		($gene_name, $cds_id) = ($refseq_info, 0);
	}

	# now do the insertion
	#
	$sql = "INSERT INTO $db_gene VALUES (0, '$items[0]', $items[1], $items[2], $cds_id, 0, 0, 0, 0, '$gene_symbol', '$gene_name')";
	$sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";

	# now need to check if the gene_name has changed
	#
	if ($prev_gene_name ne $gene_name) {
		if ($prev_gene_name ne ".") {
			# need to update the previous gene's information in the database
			#
			$sql = "UPDATE $db_gene set exon_count=$exon_count, cds_start=$cds_start, cds_end=$cds_end, cds_length=$cds_length where gene_name='$prev_gene_name'";
			$sth = $dbh->prepare($sql) or die "Query problem $!\n";
			$sth->execute() or die "Execution problem $!\n";
		} 
		
		# update all the following variable to the new gene
		# According to the CCDS website (eg, for gene NM_001077), the coordinates are in bed format
		# Therefore, we should NOT use $cds_end - $cds_start + 1;
		#
		$cds_start = $items[1];
		$cds_end = $items[2];
		$cds_length = $cds_end - $cds_start;
		$exon_count = 1;
		$prev_gene_name = $gene_name;
	} else {
		$exon_count += 1;
		$cds_length += $items[2] - $items[1];

		if ($cds_start > $items[1]) {
			$cds_start = $items[1];
		}

		if ($cds_end < $items[2]) {
			$cds_end = $items[2];
		}
	}
}
close(IN);

# now need to update cds info for the last gene
#
$sql = "UPDATE $db_gene set exon_count=$exon_count, cds_start=$cds_start, cds_end=$cds_end, cds_length=$cds_length where gene_name='$prev_gene_name'";
$sth = $dbh->prepare($sql) or die "Query problem $!\n";
$sth->execute() or die "Execution problem $!\n";

print("Finish DB dumping\n");

$dbh->disconnect();

exit;
