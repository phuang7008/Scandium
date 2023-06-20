#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

$sql = "select refseq_name, count(distinct exon_id) from Gene_RefSeq_Exon2 group by refseq_name";
$sth=$dbh->prepare($sql) or die "DB query error: $!";
$sth->execute() or die "DB execution error: $!";

my %refseq;
while ( my($refseq_name, $count) = $sth->fetchrow_array) {
	$refseq{$refseq_name} = $count;
}

$sql = "select refseq_name, count(exon_id) from Gene_RefSeq_Exon2 group by refseq_name";
$sth=$dbh->prepare($sql) or die "DB query error: $!";
$sth->execute() or die "DB execution error: $!";

my %refseq_dup;
while ( my ($refseq_name, $count) = $sth->fetchrow_array) {
	$refseq_dup{$refseq_name} = $count;
}

# now do the comparison
foreach my $ref (keys %refseq) {
	if ($refseq{$ref} != $refseq_dup{$ref}) {
		print("Duplication found in $ref\n");
	} else {
		print("Same\n");
	}
}

$dbh->disconnect();
close(IN);

exit;
