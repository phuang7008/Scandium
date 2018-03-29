#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the  annotation\n";
my $type = shift || die "Please enter the version of annotation hg38 or hg37\n";

my $db_anno = $type=~/hg38/i ? "VCRomePKv2_Annotation38" : "VCRomePKv2_Annotation37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS $db_anno");
};

# now create table again
#
$dbh->do(qq{
  CREATE TABLE $db_anno (
  `aid`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `start`  INT UNSIGNED NOT NULL,
  `end`    INT UNSIGNED NOT NULL,
  `gene_symbol`  varchar(250) NULL,
  `Synonymous`   varchar(2500) NULL,
  `prev_gene_symbol` varchar(2500) NULL,
  `annotation` varchar(4000) NOT NULL,
  PRIMARY KEY (aid),
  INDEX `CHR` (`chrom`),
  INDEX START (start),
  INDEX END (end),
  INDEX `COMP` (`chrom`, start, end),
  INDEX `SYM` (`gene_symbol`),
  INDEX `ANNO` (`annotation`)
) ENGINE=InnoDB;
});

open (IN, $file) or die "open failed for reading: $!";
while (<IN>) {
	chomp;
	next if $_!~/\w+/;

	my @items = split "\t";

	# process the 4rd (index at 3) item on the list that contains gene_symbol, refseq_name
	# 1       69090   70008   OR4F5=CCDS30547_exon_0;OR4F5=NM_001005484_exon_0;OR4F5=OTTHUMT00000003223_exon_0
	#
	my @list = split(";", $items[3]);
	my (%annotations, %gene_name_list, @all_transcripts);
	foreach my $info (@list) {
		my ($cur_gene_symbol, $refseq_info, $cur_Synonymous) = split(/\=/, $info);
		if ($cur_gene_symbol ne ".") {
			$gene_name_list{$cur_gene_symbol}++;
		}

		if ($cur_Synonymous ne ".") {
			$gene_name_list{$cur_Synonymous}++;
		}

		# push @{$annotations{$gene_symbol}}, $refseq_info;
		push @all_transcripts, $refseq_info;
	}

	my ($gene_symbol, $Synonymous);
	foreach my $g (sort {$gene_name_list{b} <=> $gene_name_list{a}} keys %gene_name_list) {
		if (!defined $gene_symbol) {
			$gene_symbol = $g;
		} else {
			if (!defined $Synonymous) {
				$Synonymous = $g;
			} else {
				$Synonymous = $Synonymous.";".$g;
			}
		}
	}

	if (!defined $gene_symbol) {
		$gene_symbol = ".";
	}

	if (!defined $Synonymous) {
		$Synonymous = ".";
	}

	# now do the insertion
	#
	my ($refseq, $ccds, $vega, $miRNA) = (".", ".", ".", ".");
	foreach	my $xp (sort @all_transcripts) {
		if ($xp=~/^N.*/) {
			if ($refseq eq ".") {
				$refseq = $xp;
			} else {
				$refseq = $refseq.";".$xp;
			}
		}

		if ($xp=~/^CCDS.*/) {
			if ($ccds eq ".") {
				$ccds = $xp;
			} else {
				$ccds = $ccds.";".$xp;
			}
		}

		if ($xp=~/^OTT.*/) {
			if ($vega eq ".") {
				$vega = $xp;
			} else {
				$vega = $vega.";".$xp;
			}
		}

		if ($xp=~/^hsa.*/) {
			if ($miRNA eq ".") {
				$miRNA = $xp;
			} else {
				$miRNA = $miRNA.";".$xp;
			}
		}
	}

	my $exon_list = $refseq."\t".$ccds."\t".$vega."\t".$miRNA."\t.";
	chomp $exon_list;

	$sql = "INSERT INTO $db_anno VALUES (0, '$items[0]', $items[1], $items[2], '$gene_symbol', '$Synonymous', '.', '$exon_list')";
	$sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";
}
close(IN);

print("Finish DB dumping\n");

$dbh->disconnect();

exit;
