#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the RefSeq annotation\n";
my $type = shift || die "Please enter the version of annotation hg38 or hg37\n";
my $flag = shift || die "Please specify if this is for gene or snp\n";

my $db_gene = $type=~/hg38/i ? "Right10K_CDS38" : "Right10K_CDS37";
my $db_snp  = $type=~/hg38/i ? "Right10K_SNP38" : "Right10K_SNP37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS $db_gene") if $flag!~/snp/i;
	$dbh->do("DROP TABLE IF EXISTS $db_snp")  if $flag=~/snp/i;
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
}) if $flag!~/snp/i;

$dbh->do(qq{
  CREATE TABLE $db_snp (
  `id`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `snp_start` INT UNSIGNED NOT NULL,
  `snp_end` INT UNSIGNED NOT NULL,
  `gene_symbol` varchar(200) NOT NULL,
  `gene_name` varchar(200) NOT NULL,
  PRIMARY KEY (id),
  INDEX `CHR` (`chrom`),
  INDEX `COMP` (`chrom`, snp_start, snp_end),
  INDEX `SYM` (`gene_symbol`),
  INDEX `NAME` (`gene_name`)
) ENGINE=InnoDB;
}) if $flag=~/snp/i;

my $prev_gene_name=".";
my $gene_name="";
my $cds_length=0;
my $exon_count=0;

open (IN, $file) or die "open failed for reading: $!";

while (<IN>) {
	chomp;
	my @items = split "\t";

	if ($items[3]=~/ENST/) {
		$items[1] -= 1;
	}

	# process the 3rd item on the list that contains the cds_id, gene name and gene symbol
	#
	my ($gene_symbol, $cds_info) = split(/\|/, $items[3]);
	my $cds_id=0;

	if ($cds_info=~/cds/i) {
		($gene_name, $cds_id) = split(/_cds_/, $cds_info);
	} else {
		($gene_name, $cds_id) = ($cds_info, 0);
	}

	if ($gene_name eq $prev_gene_name) {
		# refseq/CCDS/vega is 0 based, so don't need to + 1
		#
		$cds_length += $items[2] - $items[1] + 1;
		$exon_count += 1;
	} else {
		if ($prev_gene_name ne ".") {
			if ($flag!~/snp/i) {
				# if ($prev_gene_name=~/^NM/ || $prev_gene_name=~/^ENST/) 
				# only need to update the refseq or gencode (ensembl) gene's cds_length info for the corresponding rows
				#
				$sql = "UPDATE $db_gene SET exon_count=$exon_count, cds_length=$cds_length  WHERE gene_name='$prev_gene_name'";
				$sth = $dbh->prepare($sql) or die "Query problem $!\n";
				$sth->execute() or die "Execution problem $!\n";
			}
		}

		# now re-start the new cds_length and exon_count info here
		#
		$cds_length = $items[2] - $items[1] + 1;
		$exon_count = 1;
		$prev_gene_name = $gene_name;
	}

	# now do the insertion
	#
	if ($flag=~/snp/i) {
		# for SNP
		#
		$sql = "INSERT INTO $db_snp VALUES (0, '$items[0]', $items[1], $items[2], '$gene_symbol', '$gene_name')";
		$sth = $dbh->prepare($sql) or die "Query problem $!\n";
		$sth->execute() or die "Execution problem $!\n";
	} else {
		# for gene cds
		#
		$sql = "INSERT INTO $db_gene VALUES (0, '$items[0]', $items[1], $items[2], $cds_id, $exon_count, 0, 0, 0, '$gene_symbol', '$gene_name')";
		$sth = $dbh->prepare($sql) or die "Query problem $!\n";
		$sth->execute() or die "Execution problem $!\n";
	}
}
close(IN);

# now update the last gene cds_length info
#
if ($flag!~/snp/i) {
	$sql = "UPDATE $db_gene SET exon_count=$exon_count, cds_length=$cds_length  WHERE gene_name='$gene_name'";
	$sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";
}

print("Finish DB dumping\n");

$dbh->disconnect();

exit;
