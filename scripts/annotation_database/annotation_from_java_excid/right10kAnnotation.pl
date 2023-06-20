#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift || die "Please enter the name of the RefSeq annotation\n";
my $type = shift || die "Please enter the version of annotation hg38 or hg37\n";

my $db_anno = $type=~/hg38/i ? "Right10K_Annotation38" : "Right10K_Annotation37";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS $db_anno");
};

# now create table nameain
$dbh->do(qq{
  CREATE TABLE $db_anno (
  `aid`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `start`  INT UNSIGNED NOT NULL,
  `end`    INT UNSIGNED NOT NULL,
  `gene_symbol`  varchar(250) NULL,
  `Synonymous`   varchar(2500) NULL,
  `prev_gene_symbol` varchar(2500) NULL,
  `annotation` varchar(3000) NOT NULL,
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

    # process the 3rd item on the list that contains gene_symbol, refseq_name
    # eg, 3       123451742       123451949       MYLK|NM_053027_cds_22;MYLK|NM_053025_cds_23
    #
    my @list = split(";", $items[3]);
    my %annotations;

    foreach my $info (@list) {
        my ($gene_symbol, $refseq_info) = split(/\|/, $info);
        push @{$annotations{$gene_symbol}}, $refseq_info;
    }

    # now do the insertion
    #
    foreach my $g_symbol (keys %annotations) {
		# my $cds_list = join(";", sort @{$annotations{$g_symbol}});
		#
		my $cds_list="0";
		my $refseq_list  = ".";
		my $gencode_list = ".";
		my $other_list =".";

		foreach my $name (sort @{$annotations{$g_symbol}}) {
			if ($cds_list eq "0") {
				$gencode_list = $name if ($name=~/ENST/i);
				$refseq_list  = $name if ($name=~/NM|NR/i);
				$other_list   = $name if ($name=~/SNP|Pseudo/i);
			} else {
				if ($name=~/ENST/i) {
					if (length($gencode_list) > 1) {
						$gencode_list .= ";$name";
					} else {
						$gencode_list = $name;
					}
				}

				if ($name=~/NM|NR/i) {	
					if (length($refseq_list) > 1) {
						$refseq_list .= ";$name";
					} else {
						$refseq_list = $name;
					}
				}

				if ($name=~/SNP|Pseudo/i) {
					if (length($other_list) > 1) {
						$other_list .= ";$name";
					} else {
						$other_list = $name;
					}
				}
			}
		}

		$cds_list = $refseq_list."\t.\t".$gencode_list."\t.".$other_list."\n";

        $sql = "INSERT INTO $db_anno VALUES (0, '$items[0]', $items[1], $items[2], '$g_symbol', '.', '.', '$cds_list')";
        $sth = $dbh->prepare($sql) or die "Query problem $!\n";
        $sth->execute() or die "Execution problem $!\n";
    }
}
close(IN);

print("Finish DB dumping\n");

$dbh->disconnect();

exit;
