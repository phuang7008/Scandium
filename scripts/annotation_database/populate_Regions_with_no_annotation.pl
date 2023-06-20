#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;
use DBI;

my $file = shift or die "Please enter the file name that contains all inter-genic regions\n";

# connect to the database
my $dbh = DBI->connect('DBI:mysql:GeneAnnotations:sug-esxa-db1', 'phuang', 'phuang') or die "DB connection failed: $!";

my ($sql, $sth);

# drop the table first
eval {
	$dbh->do("DROP TABLE IF EXISTS Intergenic_Regions");
};

# now create table again
$dbh->do(qq{
CREATE TABLE Intergenic_Regions (
  `no_id`  INT UNSIGNED NOT NULL AUTO_INCREMENT,
  `chrom`  varchar(50) NOT NULL,
  `start`  INT UNSIGNED NOT NULL,
  `end`    INT UNSIGNED NOT NULL,
  PRIMARY KEY (no_id),
  INDEX `CHR` (`chrom`),
  INDEX START (start),
  INDEX END (end),
  INDEX `COMP` (`chrom`, start, end)
) ENGINE=InnoDB;
});

open (IN, $file) or die "file open failed: $!\n";

while (<IN>) {
	chomp;
	my @items = split "\t";

	$sql = "INSERT INTO Intergenic_Regions VALUES(0, '$items[0]', $items[1], $items[2])";
	$sth = $dbh->prepare($sql) or die "DB query error: $!";
	$sth->execute() or die "DB execution error: $!";
}

$dbh->disconnect();
close(IN);

exit;
