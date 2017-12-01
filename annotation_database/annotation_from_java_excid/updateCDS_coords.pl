#!/hgsc_software/perl/perl-5.18.2/bin/perl

# This script is used to update CDS coordinates for all gene/exon/transcript annotation databases
# Before running this script, we need to obtain the CDS coordinates from all the gene annotation sources
# such as RefSeq, CCDS, Gencode, VEGA etc.
# The miRNA info doesn't need to be changed as they only contain one exon, and being handled already!
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the file name of the sorted gene annotation\n";
my $db = shift || die "Please enter the name of database to be modified\n";

#my $db_eMerge      = $type=~/hg38/i ? "eMerge_CDS38" : "eMerge_CDS37";
#my $db_right10K    = $type=~/hg38/i ? "Right10K_CDS38" : "Right10K_CDS37";
#my $db_VCRome_PKv2 = $type=~/hg38/i ? "Gene_RefSeq_CDS38" : "Gene_RefSeq_CDS37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);
my ($prev_cds_start, $prev_cds_end);

my (%refseqCoords, %ccdsCoords, %vegaCoords, %gencodeCoords);
my $prev_gene_name = '.';

open (IN, $file) or die "open failed for reading: $!";

while (<IN>) {
	chomp;
	my @items = split "\t";

	# process the 3rd item on the list that contains the 
	# gene name, cds_id, exon_count, gene symbol, cds_start, cds_end, cds_length
	# eg, CCDS2_exon_1_13=SAMD11=861321=879533=2046
	# eg, NM_000016_exon_1_12=ACADM=76190472=76228448=1266
	#
	my ($cds_info, $gene_symbol, $cds_start, $cds_end, $cds_length) = split(/=/, $items[3]);
	my ($gene_name, , $cds_id) = split(/_exon_/, $cds_info);
	if ($gene_name=~/(ENST\d+)\.\d+/) {
		$gene_name=$1;
	}

	if ($gene_name eq $prev_gene_name) {
		next;
	} else {
		# start a new gene cds, do update...
		#
		if ($prev_gene_name ne ".") {
			if ($db eq "eMerge_CDS37" || $db eq "Right10K_CDS37") {
				$sql = "UPDATE $db SET cds_start=$prev_cds_start, cds_end=$prev_cds_end  WHERE refseq_name='$prev_gene_name'";
			} else {
				$sql = "UPDATE $db SET cds_start=$prev_cds_start, cds_end=$prev_cds_end  WHERE gene_name='$prev_gene_name'";
			}
			$sth = $dbh->prepare($sql) or die "Query problem $!\n";
			$sth->execute() or die "Execution problem $!\n";
		}

		$prev_gene_name = $gene_name;
		$prev_cds_start = $cds_start;
		$prev_cds_end   = $cds_end;
	}
}
close(IN);

# process the last one
#
if ($db eq "eMerge_CDS37" || $db eq "Right10K_CDS37") {
	$sql = "UPDATE $db SET cds_start=$prev_cds_start, cds_end=$prev_cds_end  WHERE refseq_name='$prev_gene_name'";
} else {
	$sql = "UPDATE $db SET cds_start=$prev_cds_start, cds_end=$prev_cds_end  WHERE gene_name='$prev_gene_name'";
}

$sth = $dbh->prepare($sql) or die "Query problem $!\n";
$sth->execute() or die "Execution problem $!\n";

print("Finish DB Updating\n");

$dbh->disconnect();

exit;
