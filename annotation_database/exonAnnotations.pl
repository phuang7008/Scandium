#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the file that contains all the intronic regions\n";
my $type = shift || die "Please specify the version of gene annotation hg19 or hg38\n";
my $database = $type eq "hg38" ? "Exon_Regions38" : "Exon_Regions37";
my $hgnc = $type eq "hg38" ? "HGNC38" : "HGNC37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS $database");
};

# now create table again
$dbh->do(qq{
CREATE TABLE $database (
  `ex_id`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `start`  INT UNSIGNED NOT NULL,
  `end`    INT UNSIGNED NOT NULL,
  `gene_symbol`  varchar(250) NULL,
  `Synonymous`   varchar(2500) NULL,
  `prev_gene_symbol` varchar(2500) NULL,
  `annotation` varchar(3000) NOT NULL,
  PRIMARY KEY (ex_id),
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
	my @items = split "\t";

	# process the 4th entry that contains the gene info
	my @info = split(/;/, $items[3]);

	my %genes;
	my %prev_genes;
	my %exons;

	foreach my $combo (@info) {
		my ($name, $name2) = split(/=/, $combo);
		$genes{$name2}++ if (defined $name2 and $name2 ne "");
		$exons{$name}++;
	}

	foreach my $sn (keys %exons) {
		$sn=~s/\_\d+\_\d+$//;

		my $sql;

        if ($sn=~/^N.*/ || $sn=~/^X.*/) {
            $sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (refseq_accession IS NOT NULL AND find_in_set('$sn', replace(refseq_accession, '|', ',')))";
        }

		if ($type eq "hg38") {
			if ($sn=~/^uc.*/) {
				$sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (ucsc_id IS NOT NULL AND find_in_set('$sn', replace(ucsc_id, '|', ',')))";
			}

			if ($sn=~/^EN.*/) {
				$sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (ensembl_gene_id IS NOT NULL AND find_in_set('$sn', replace(ensembl_gene_id, '|', ',')))";
			}
		} else {
			if ($sn=~/^OTT.*/) {
				$sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (vega_id IS NOT NULL AND find_in_set('$sn', replace(vega_id, '|', ',')))";
			}
        }

        if ($sn=~/^CCDS.*/) {
            $sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (ccds_id IS NOT NULL AND find_in_set('$sn', replace(ccds_id, '|', ',')))";
        }

		if ($sn=~/^hsa.*/) {
			$sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (mirbase IS NOT NULL AND find_in_set('$sn', replace(mirbase, '|', ',')))";
		}

		#my $sql = "SELECT symbol, prev_symbol FROM HGNC WHERE (symbol='$name2') OR (symbol='$name') OR (ccds_id IS NOT NULL AND find_in_set('$name', replace(ccds_id, '|', ','))) OR (refseq_accession IS NOT NULL AND find_in_set('$name', replace(refseq_accession, '|', ','))) OR (vega_id IS NOT NULL AND find_in_set('$name', replace(vega_id, '|', ',')))";
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

	my ($refseq, $ccds, $vega, $gencode, $miRNA) = (".", ".", ".", ".", ".");

	foreach my $ex (sort keys %exons) {
		next if ($ex eq "." || $ex eq "");
		$ex=~s/\_\d+$//;

		if ($ex=~/^N.*/ || $ex=~/^X.*/) {
			$refseq = $refseq eq "." ? $ex : "$refseq;$ex";
		}

		if ($ex=~/^CCDS.*/) {
			$ccds = $ccds eq "." ? $ex : "$ccds;$ex";
		}

		if ($ex=~/^OTT.*/) {
			$vega = $vega eq "." ? $ex : "$vega;$ex";
		}

		if ($ex=~/^uc.*/) {
			$gencode = $gencode eq "." ? $ex : "$gencode;$ex";
		}

		if ($ex=~/^EN.*/) {
            $gencode = $gencode eq "." ? $ex : "$gencode;$ex";
        }

		if ($ex=~/^hsa.*/) {
			$miRNA = $miRNA eq "." ? $ex : "$miRNA;$ex";
		}
	}

	my $annotations = $type eq "hg38" ? "$refseq\t$ccds\t$gencode\t$miRNA" : "$refseq\t$ccds\t$vega\t$miRNA";

	my $sql = "INSERT INTO $database VALUES (0, '$items[0]', $items[1], $items[2], '$gene', '$Synonymous', '$prev_gene', '$annotations')";
	my $sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";
}

$dbh->disconnect();
close(IN);

exit;
