#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use FindBin;
use lib "$FindBin::Bin";
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the file that contains all the exon regions\n";
my $type = shift || die "Please specify the version of gene annotation hg19 or hg38\n";
my $hgnc = $type eq "hg38" ? "HGNC38" : "HGNC37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang468', 'phuang468') or die "DB connection failed: $!";

my ($sql, $sth);

# Here is the detailed item list
# 0: transcript name; 1: chromosome id; 2: strand orientation; 3: transcript start; 4 transcript end; 
# 5: cds start; 6: cds end; 7: exon count; 8: exon starts; 9: exon ends; 10: gene symbol
#
my %refseq_names;

open (IN, $file) or die "open failed for reading: $!";
while (<IN>) {
	next if /chrom/i;
	chomp;
	my @items = split "\t";

	my @starts = split ",", $items[8];
	my @ends   = split ",", $items[9];

	# remove 'chr' in front of chromosome id
	#$items[1]=~s/chr//i;

	if ($items[0]=~/^CCDS.*/) {
		$items[0]=~s/(CCDS.*)\.\d+$/$1/;
		$sql = "SELECT symbol FROM $hgnc WHERE (ccds_id IS NOT NULL AND find_in_set('$items[0]', replace(ccds_id, '|', ',')))";
		$sth = $dbh->prepare($sql) or die "DB query error: $!";
		$sth->execute() or die "DB execution error: $!";

		my ($gene_symbol) = $sth->fetchrow_array;
		$items[10] = $gene_symbol;
		$sth->finish();
	}

	if (! defined $items[10]) { $items[10] = "."; }

	# Some refseq names are the same even if they point to two different transcripts at different genomic location
	# Here we need to separate them 
	#
	if ($items[0]=~/^N.*/) {
		if (defined $refseq_names{$items[0]}) {
			my $tmp_name = $items[0]."-$refseq_names{$items[0]}";
			$refseq_names{$items[0]}++;
			$items[0] = $tmp_name;
		} else {
			$refseq_names{$items[0]}++;
		}
	}

	my $cds_start = $items[5];
	my $cds_end   = $items[6];

	# we need to loop through the exon regions twice. The first time to calculate the cds size    
	# while the second time will be to print the combined results out to a file                   
	# for 'NR' genes, the cds start = cds end, so we will use transcript start and end!           
	#
	if ($cds_start == $cds_end) {
		$cds_start = $items[3];                                                                 
		$cds_end   = $items[4];
	}

	my $cds_length = 0;

	foreach my $idx (0..$#starts) {
		# exon_start ========= exon end                                                           
		#                     cds_start ---------- cds_end                                        
		#                                                                                         
		next if ($ends[$idx] < $cds_start);                                                          

		#               exon_start ========== exon_end                                            
		# cds_start -------- cds_end                                                              
		#                                                                                         
		next if ($starts[$idx] > $cds_end); 

		# exon_start ======= exon_end                                                             
		#     cds_start -------- cds_end                                                          
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

	# second loop will output the combined results
	#
	foreach my $idx (0..$#starts) {

		my $exon_id = $idx;
		if ($items[2] eq "-") {
			$exon_id = $items[7] - $idx - 1;
		}

		#print "$items[1]\t$starts[$idx]\t$ends[$idx]\t$items[0]_exon_$exon_id\_$items[7]=$items[10]\t$items[10]\n";
		#print "$items[1]\t$starts[$idx]\t$ends[$idx]\t$items[0]_exon_$exon_id\_$items[7]=$items[10]=$cds_start=$cds_end=$cds_length\t$items[10]\n";
		print "$items[1]\t$starts[$idx]\t$ends[$idx]\t$items[0]|exon_$exon_id|gene\_$items[7]=$items[10]=$cds_start=$cds_end=$cds_length\t$items[10]\n";
	}
}

$dbh->disconnect();
close(IN);

exit;
