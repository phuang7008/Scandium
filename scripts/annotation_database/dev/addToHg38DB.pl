#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the RefSeq annotation\n";
my $type = shift || die "Please enter the version of annotation hg38 or hg37\n";

my $database = $type=~/hg38/i ? "Gene_RefSeq_Dynamic_CDS_Coords38" : "Gene_RefSeq_Dynamic_CDS_Coords37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# Here is the sample 
# name   chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        name2
# NR_046430.1     chr4    +       75556047        75565893        75565893        75565893        4       75556047,75563209,75564113,75564387,    75556149,75563345,75564252,75565893,    C4orf26
#
open (IN, $file) or die "open failed for reading: $!";
my $header=<IN>;

while (<IN>) {
	chomp;
	my @items  = split "\t";
	my @starts = split(/\,/, $items[8]);
	my @ends   = split(/\,/, $items[9]);
	my $cds_start = $items[5];
	my $cds_end   = $items[6];
	my @names = split(/\./, $items[0]);
	$items[0]=~s/\.\d+$//;
	$items[1]=~s/chr//i;
	my $cds_length = 0;

	foreach my $idx (0..$#starts) {
		# need to cut off non-coding regions
		# exon_start ========= exon end
		#							cds_start ---------- cds_end  
		#
		next if ($ends[$idx] < $cds_start);

		#               exon_start ========== exon_end 
		# cds_start -------- cds_end
		#
		next if ($starts[$idx] > $cds_end);

		# exon_start ======= exon_end 
		#      cds_start --------- cds_end 
		#
		if ( $starts[$idx] < $cds_start and $ends[$idx] < $cds_end) {
			$cds_length += $ends[$idx] - $cds_start;
			next;
		}

		#      exon_start ======== exon_end
		# cds_start --------- cds_end 
		#
		if ($cds_start < $starts[$idx] and $cds_end < $ends[$idx]) {
			$cds_length += $cds_end - $starts[$idx];
			next;
		}

		$cds_length += $ends[$idx] - $starts[$idx];
	}

	 # second loop will 
	 #
	 foreach my $idx (0..$#starts) {
		 next if ($ends[$idx] < $cds_start);
		 next if ($starts[$idx] > $cds_end);
		 $starts[$idx] = $cds_start if ($starts[$idx] < $cds_start);
		 $ends[$idx]   = $cds_end if ($ends[$idx] > $cds_end);

		 my $exon_id = $idx;
		 if ($items[2] eq "-") {
			 $exon_id = $items[7] - $idx - 1;
		 }

		my ($sql, $sth);
		if ($ends[$idx] > $starts[$idx]) {
			# now do the insertion
			#
			$sql = "INSERT INTO $database VALUES (0, '$items[1]', $starts[$idx], $ends[$idx], $exon_id, $items[7], $cds_start, $cds_end, $cds_length, '$items[10]', '$items[0]')";
			$sth = $dbh->prepare($sql) or die "Query problem $!\n";
			$sth->execute() or die "Execution problem $!\n";
			#print("$sql\n");
		}
	}
}

print("Finish adding database to DB\n");

$dbh->disconnect();
close(IN);

exit;
