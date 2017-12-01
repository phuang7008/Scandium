#!/hgsc_software/perl/perl-5.18.2/bin/perl

# This script is used to update Database with cds_start and cds_end entries with 0s.
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $db = shift || die "Please enter the name of database to be modified\n";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

#my $sql = "select distinct(gene_name) from Gene_RefSeq_CDS37 where cds_start=0";
my $sql = "select distinct(gene_name) from $db where cds_start=0";
my $sth = $dbh->prepare($sql) or die "Query problem $!\n";                                           
$sth->execute() or die "Execution problem $!\n";

my %gene_names;
while ( my ($tmp_name) = $sth->fetchrow_array) {
	$gene_names{$tmp_name}++;
}

my ($prev_cds_start, $prev_cds_end) = 0;

foreach my $name (keys %gene_names) {
	$sql = "select cds_target_start, cds_target_end from $db where gene_name='$name'";
	$sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";

	while ( my($cds_start, $cds_end) = $sth->fetchrow_array) {
		if ($cds_start < $prev_cds_start or $prev_cds_start == 0) {
			$prev_cds_start = $cds_start;
		}

		if ($cds_end > $prev_cds_end or $prev_cds_end == 0) {
			$prev_cds_end = $cds_end;
		}
	}

	# update here
	#
	$sql = "UPDATE $db SET cds_start=$prev_cds_start, cds_end=$prev_cds_end  WHERE gene_name='$name'";
	$sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";

	# reset 
	#
	$prev_cds_start = 0;
	$prev_cds_end = 0;
}

print("Finish DB Updating\n");

$dbh->disconnect();

exit;
