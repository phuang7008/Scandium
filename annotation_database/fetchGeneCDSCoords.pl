#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the RefSeq annotation\n";
my $type = shift || die "Please enter the version of annotation hg38 or hg37\n";
my $project_name = shift || die "Please enter the name of the project, such as eMerge, right10K or VCRome_PKv2\n";

my $database = $type=~/hg38/i ? "Gene_".$project_name."_CDS_Coords38" : "Gene_".$project_name."_CDS_Coords37";
my $hgnc     = $type=~/hg38/i ? "HGNC38" : "HGNC37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS $database");
};

# now create table again
# cds_target_start ==> cds exon start (so it is exon based)
# cds_start => cds transcript start (so it is the entire CDS)
#
$dbh->do(qq{
  CREATE TABLE $database (
  `id`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `cds_target_start` INT UNSIGNED NOT NULL,
  `cds_target_end` INT UNSIGNED NOT NULL,
  `exon_id` INT UNSIGNED NOT NULL,
  `exon_count` INT UNSIGNED NOT NULL,
  `cds_start`  INT UNSIGNED NOT NULL,
  `cds_end`    INT UNSIGNED NOT NULL,
  `cds_length` INT SIGNED NOT NULL,
  `gene_symbol` varchar(200) NOT NULL,
  `transcript_name`   varchar(200) NOT NULL,
  PRIMARY KEY (id),
  INDEX `CHR` (`chrom`),
  INDEX START (cds_start),
  INDEX END (cds_end),
  INDEX `COMP` (`chrom`, cds_target_start, cds_target_end),
  INDEX `SYM`  (`gene_symbol`),
  INDEX `NAME` (`transcript_name`)
) ENGINE=InnoDB;
});

open (IN, $file) or die "open failed for reading: $!";
#my $header=<IN>;

while (<IN>) {
	chomp;
	my @items = split "\t";

	# process the 4rd item (index 3) on the list that contains the exon_id, refseq name and gene symbol
	# 1       868070  868675  NR_027055_2_3=FAM41C=868070=876802=1706
	# 1       30365   30503   ENST00000607096.1_0_1=MIR1302-11=30365=30503=138
	#
	my @info = split(/;/, $items[3]);

	foreach my $region (@info) {
		my ($ref_exon_id, $gene_symbol, $cds_start, $cds_end, $cds_length) = split(/=/, $region);
		my ($refseq_name, $exon_id, $exon_count);
		#if ($ref_exon_id=~/(^N.*)_(\d+)_(\d+)$/) {
	   	if ($ref_exon_id=~/(^N.*)_(\d+)_(\d+)$/) {
			($refseq_name, $exon_id, $exon_count) = ($1, $2, $3);
		} elsif ($ref_exon_id=~/(^X.*)_(\d+)_(\d+)$/) {
			($refseq_name, $exon_id, $exon_count) = ($1, $2, $3);
		} elsif ($ref_exon_id=~/(^E.*)_(\d+)_(\d+)$/) {
			($refseq_name, $exon_id, $exon_count) = ($1, $2, $3);
			$refseq_name=~s/\.\d+//;
		}
		
		my ($sql, $sth);
		if (!defined $gene_symbol) {
			$sql = "SELECT symbol FROM $hgnc WHERE (refseq_accession IS NOT NULL AND find_in_set('$refseq_name', refseq_accession))";
			$sth = $dbh->prepare($sql) or die "DB query error: $!";
			$sth->execute() or die "DB execution error: $!";

			($gene_symbol) = $sth->fetchrow_array if (!defined $gene_symbol);
		}

		# we only need CDS information. So need to check to exclude UTRs from CDS
		# skip it is current exon is all UTRs
		#
		next if ($items[2] < $cds_start);
		
		$items[1] = $cds_start if ($items[1] < $cds_start);
		$items[2] = $cds_end   if ($items[2] > $cds_end);

		# now do the insertion
		#
		$sql = "INSERT INTO $database VALUES (0, '$items[0]', $items[1], $items[2], $exon_id, $exon_count, $cds_start, $cds_end, $cds_length, '$gene_symbol', '$refseq_name')";
		$sth = $dbh->prepare($sql) or die "Query problem $!\n";
		$sth->execute() or die "Execution problem $!\n";
	}
}

print("Finish DB dumping\n");

# now I need to update the exon_count information
#$sql = "select refseq_name, count(exon_id) from $database group by refseq_name";
#$sth=$dbh->prepare($sql) or die "DB query error: $!";
#$sth->execute() or die "DB execution error: $!";

#my %refseq;
#while ( my($refseq_name, $count) = $sth->fetchrow_array) {
#    $refseq{$refseq_name} = $count;
#}

#foreach my $refseq_name (keys %refseq) {
    # now update the table
#    $sql = "UPDATE $database SET exon_count=$refseq{$refseq_name} WHERE refseq_name='$refseq_name'";
#    $sth = $dbh->prepare($sql) or die "Query problem $!\n";
#    $sth->execute() or die "Execution problem $!\n";
#}

$dbh->disconnect();
close(IN);

exit;
