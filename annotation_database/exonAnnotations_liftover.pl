#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the file that contains all the intronic regions\n";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS Capture_Regions38");
};

# now create table again
$dbh->do(qq{
CREATE TABLE Capture_Regions38 (
  `ex_id`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `start`  INT UNSIGNED NOT NULL,
  `end`    INT UNSIGNED NOT NULL,
  `exon_start`  INT UNSIGNED NOT NULL,
  `exon_end`    INT UNSIGNED NOT NULL,
  `capture_start`  INT UNSIGNED NOT NULL,
  `capture_end`    INT UNSIGNED NOT NULL,
  `gene_symbol`  varchar(250) NULL,
  `Synonymous`   varchar(2500) NULL,
  `prev_gene_symbol` varchar(2500) NULL,
  `annotation` varchar(3000) NOT NULL,
  PRIMARY KEY (ex_id),
  INDEX `CHR` (`chrom`),
  INDEX START (start),
  INDEX END (end),
  INDEX EXON_START (exon_start),
  INDEX EXON_END (exon_end),
  INDEX CAP_START (capture_start),
  INDEX CAP_END (capture_end),
  INDEX `COMP` (`chrom`, start, end),
  INDEX `SYM` (`gene_symbol`),
  INDEX `ANNO` (`annotation`)
) ENGINE=InnoDB;
});

open (IN, $file) or die "open failed for reading: $!";

while (<IN>) {
	chomp;
	my @items = split "\t";

	# process the 4th entry that contains the gene info
	my @info = split(/;/, $items[3]);
	#my @info = split(/,/, $items[4]);

	my %genes;
	my %prev_genes;
	my %exons;
	my %exon_starts;
	my %exon_ends;

	foreach my $combo (@info) {
		my ($name, $name2, $exon_start, $exon_end) = split(/=/, $combo);
		$genes{$name2}++ if (defined $name2 and $name2 ne "");
		$exons{$name}++;
		$exons{$name}++;
		$exon_starts{$exon_start}++;
		$exon_ends{$exon_end}++;
	}

	my $ex_starts = join("|", keys %exon_starts);
	my $ex_ends   = join("|", keys %exon_ends);

	foreach my $sn (keys %exons) {
		$sn=~s/\_\d+$//;

		my $sql;

        if ($sn=~/^N.*/) {
            $sql = "SELECT symbol, prev_symbol FROM HGNC37 WHERE (refseq_accession IS NOT NULL AND find_in_set('$sn', replace(refseq_accession, '|', ',')))";
        }

        if ($sn=~/^EN.*/) {
            $sql = "SELECT symbol, prev_symbol FROM HGNC37 WHERE (ensembl_gene_id IS NOT NULL AND find_in_set('$sn', replace(ensembl_gene_id, '|', ',')))";
        }

        if ($sn=~/^CCDS.*/) {
            $sql = "SELECT symbol, prev_symbol FROM HGNC37 WHERE (ccds_id IS NOT NULL AND find_in_set('$sn', replace(ccds_id, '|', ',')))";
        }

		if ($sn=~/^hsa.*/) {
			$sql = "SELECT symbol, prev_symbol FROM HGNC37 WHERE (mirbase IS NOT NULL AND find_in_set('$sn', replace(mirbase, '|', ',')))";
		}

		my $sth = $dbh->prepare($sql) or die "DB query error: $!";
		$sth->execute() or die "DB execution error: $!";

		while ( my($gene_symbol, $prev_gene_symbol) = $sth->fetchrow_array) {
			$genes{$gene_symbol}++ if (defined $gene_symbol && $gene_symbol ne "");
			$prev_genes{$prev_gene_symbol}++ if (defined $prev_gene_symbol && $prev_gene_symbol ne "");
		}
	}

	# now produce compound gene and prev_gene info
	my $gene=".";
	my $Synonymous=".";
	foreach my $gn (sort keys %genes) {
		next if ($gn eq "." || $gn eq "");

		if ($gene eq ".") {
			$gene = $gn;
		} else {
			if ($Synonymous eq ".") {
				$Synonymous = $gn;
			} else {
				$Synonymous = $Synonymous.";".$gn;
			}
		}
	}

	my $prev_gene=".";
	foreach my $pg (sort keys %prev_genes) {
		next if ($pg eq "." || $pg eq "");

		if ($prev_gene eq ".") {
			$prev_gene = $pg;
		} else {
			$prev_gene = $prev_gene.";".$pg;
		}
	}

	my ($refseq, $ccds, $gencode, $miRNA) = (".", ".", ".", ".");
	if ($items[1] == 17369) {
		print("inside\n");
	}

	foreach my $ex (sort keys %exons) {
		next if ($ex eq "." || $ex eq "");

		if ($ex=~/^N.*/) {
			$refseq = $refseq eq "." ? $ex : "$refseq;$ex";
		}

		if ($ex=~/^CCDS.*/) {
			$ccds = $ccds eq "." ? $ex : "$ccds;$ex";
		}

		if ($ex=~/^EN.*/) {
			$gencode = $gencode eq "." ? $ex : "$gencode;$ex";
		}

		if ($ex=~/^hsa.*/) {
			$miRNA = $miRNA eq "." ? $ex : "$miRNA;$ex";
		}
	}

	my $annotations = "$refseq\t$ccds\t$gencode\t$miRNA";

	my $sql = "INSERT INTO Capture_Regions38 VALUES (0, '$items[0]', $items[1], $items[2], $ex_starts, $ex_ends, $items[5], $items[6], '$gene', '$Synonymous', '$prev_gene', '$annotations')";
	my $sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";
}

$dbh->disconnect();
close(IN);

exit;
