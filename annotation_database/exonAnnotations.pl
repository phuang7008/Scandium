#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the exon/gene annotation file\n";
my $type = shift || die "Please specify the version of gene annotation hg19 or hg38\n";
my $database = $type=~/hg38/i ? "Gene_Annotations38" : "Gene_Annotations37";
my $hgnc     = $type=~/hg38/i ? "HGNC38" : "HGNC37";

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
	my @items = split "\t";

	# each line looks like the following:
	# 1		11868	11873	NR_148357_exon_0_3=LOC102725121=11868=14362=1610;ENST00000456328.2_exon_0_3=DDX11L1=11868=14409=1657
	# 1		11868	11873	NR_148357|exon_0|gene_3=LOC102725121=11868=14362=1610;ENST00000456328.2|exon_0|gene_3=DDX11L1=11868=14409=1657

	# process the 4th entry that contains the gene info
	# NR_148357|exon_0|gene_3=LOC102725121=11868=14362=1610;ENST00000456328.2|exon_0|gene_3=DDX11L1=11868=14409=1657
	#
	my @info = split(/;/, $items[3]);

	my %genes;
    my %main_genes;
	my %prev_genes;
	my %exons;
    my %seen;
    #my ($refseq_gname, $ucsc_gname, $ensembl_gname, $ccds_gname, $miRNA_gname)=(".",".",".",".",".");

	foreach my $combo (sort { $b cmp $a } @info) {
        if ($seen{$combo}) { next; }

        # split the following pattern: 
        # NR_148357|exon_0|gene_3=LOC102725121=11868=14362=1610
        #
        my ($name, $gname) = split(/=/, $combo);
        $genes{$gname}++ if (defined $gname and $gname ne "");
		$exons{$name}++;
        $seen{$combo}++;
	}

	# handle exon part which is NR_148357|exon_0|gene_3
	#
	foreach my $ex (keys %exons) {
		#$ex=~s/\_exon_\d+\_\d+$//;
		my @ex_info = split(/\|/, $ex);

		my $sql;
        my ($refseq_flag, $ucsc_flag, $ensembl_flag, $ccds_flag, $miRNA_flag) = (0,0,0,0,0);

		# handle refseq gene symbols
		#
        if ($ex_info[0]=~/^N.*/ || $ex_info[0]=~/^X.*/) {
            $sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (refseq_accession IS NOT NULL AND find_in_set('$ex_info[0]', replace(refseq_accession, '|', ',')))";
            $refseq_flag = 1;
        }

		# handle UCSC gene symbols
		#
		if ($ex_info[0]=~/^uc.*/) {
			$sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (ucsc_id IS NOT NULL AND find_in_set('$ex_info[0]', replace(ucsc_id, '|', ',')))";
            $ucsc_flag=1;
		}

		# handle ENSEMBL gene symbols
		#
		if ($ex_info[0]=~/^EN.*/) {
			$sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (ensembl_gene_id IS NOT NULL AND find_in_set('$ex_info[0]', replace(ensembl_gene_id, '|', ',')))";
            $ensembl_flag=1;
		}

		# handle CCDS gene symbols
		#
        if ($ex_info[0]=~/^CCDS.*/) {
            $sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (ccds_id IS NOT NULL AND find_in_set('$ex_info[0]', replace(ccds_id, '|', ',')))";
            $ccds_flag=1;
        }

		# handle miRNA information
		#
		if ($ex_info[0]=~/^hsa.*/) {
			$sql = "SELECT symbol, prev_symbol FROM $hgnc WHERE (mirbase IS NOT NULL AND find_in_set('$ex_info[0]', replace(mirbase, '|', ',')))";
            $miRNA_flag=1;
		}

		my $sth = $dbh->prepare($sql) or die "DB query error: $!";
		$sth->execute() or die "DB execution error: $!";

		while ( my($gene_symbol, $prev_gene_symbol) = $sth->fetchrow_array) {
            if (defined $gene_symbol && $gene_symbol ne "") {
                if ($refseq_flag == 1)  { $main_genes{$gene_symbol}{"refseq"} = "$ex_info[0]|$ex_info[1]|gene"; }
                if ($ensembl_flag == 1) { $main_genes{$gene_symbol}{"ensembl"}= "$ex_info[0]|$ex_info[1]|gene"; }
                if ($ucsc_flag == 1)    { $main_genes{$gene_symbol}{"ucsc"}  = "$ex_info[0]|$ex_info[1]|gene"; }
                if ($ccds_flag == 1)    { $main_genes{$gene_symbol}{"ccds"}  = "$ex_info[0]|$ex_info[1]|gene"; }
                if ($miRNA_flag == 1)   { $main_genes{$gene_symbol}{"miRNA"} = "$ex_info[0]|$ex_info[1]|gene"; }
                
            }
			$prev_genes{$prev_gene_symbol}++ if (defined $prev_gene_symbol && $prev_gene_symbol ne "");
		}
	}

	# now produce compound gene and prev_gene info
	#
	my $gene=".";
	my $Synonymous=".";
	my ($refseq, $ccds, $gencode, $miRNA) = (".", ".", ".", ".", ".");
    my %seen_transcripts;

    foreach my $gn (keys %main_genes) {
        if (defined $main_genes{$gn}{"refseq"}) {
            if ($gene eq ".") { 
                $gene = $gn; 
                $refseq = $main_genes{$gn}{"refseq"};
                $seen_transcripts{$refseq}++;
            }
        } elsif (defined $main_genes{$gn}{"ensembl"}) {
            if ($gene eq ".") {
                $gene = $gn;
                $gencode = $main_genes{$gn}{"ensembl"};
                $seen_transcripts{$gencode}++;
            }
        } elsif (defined $main_genes{$gn}{"ucsc"}) {
            if ($gene eq ".") {
                $gene = $gn;
            }
        } elsif (defined $main_genes{$gn}{"ccds"}) {
            if ($gene eq ".") { 
                $gene = $gn; 
                $ccds = $main_genes{$gn}{"ccds"};
                $seen_transcripts{$ccds}++;
            }
        } elsif (defined $main_genes{$gn}{"miRNA"}) {
            if ($gene eq ".") { 
                $gene = $gn;
                $miRNA = $main_genes{$gn}{"miRNA"};
                $seen_transcripts{$miRNA}++;
            }
        }
    }

	foreach my $gn (sort keys %genes) {
		next if ($gn eq "." || $gn eq "");
        next if ($gn eq $gene);

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


	foreach my $ex (sort keys %exons) {
		next if ($ex eq "." || $ex eq "");
		$ex=~s/\_\d+$//;
        next if $seen_transcripts{$ex};
        $seen_transcripts{$ex}++;

		if ($ex=~/^N.*/ || $ex=~/^X.*/) {
            next if $refseq eq $ex;
			$refseq = $refseq eq "." ? $ex : "$refseq;$ex";
		}

		if ($ex=~/^CCDS.*/) {
            next if $ccds eq $ex;
			$ccds = $ccds eq "." ? $ex : "$ccds;$ex";
		}

        #if ($ex=~/^uc.*/) {
        #    next if $gencode eq $ex;
        #	$gencode = $gencode eq "." ? $ex : "$gencode;$ex";
        #}

		if ($ex=~/^EN.*/) {
            next if $gencode eq $ex;
            $gencode = $gencode eq "." ? $ex : "$gencode;$ex";
        }

		if ($ex=~/^hsa.*/) {
            next if $miRNA eq $ex;
			$miRNA = $miRNA eq "." ? $ex : "$miRNA;$ex";
		}
	}

	# the last column is for SNPs and Pseudo-genes
	#
	my $annotations = "$refseq\t$ccds\t$gencode\t$miRNA\t.";

	my $sql = "INSERT INTO $database VALUES (0, '$items[0]', $items[1], $items[2], '$gene', '$Synonymous', '$prev_gene', '$annotations')";
	my $sth = $dbh->prepare($sql) or die "Query problem $!\n";
	$sth->execute() or die "Execution problem $!\n";
}

$dbh->disconnect();
close(IN);

exit;
