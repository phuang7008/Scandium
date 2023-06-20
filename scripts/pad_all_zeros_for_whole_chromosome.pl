#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;

# this script is used to add the 0 coverage regions to the existing bed file
# 
my $in_file  = shift || die "Please specify the input bed file\n";
my $out_file = shift || die "Please specify the output file name\n";
my $chrom_id = shift || die "Please specify the chromosome id\n";
my $db_version = shift || die "Please specify the db version (such as hg19 or hg38)\n";

my %chrom_length;
if ($db_version eq "hg38") {
	%chrom_length = ("1"=>248956422, "2"=>242193529, "3"=>198295559, "4"=>190214555, "5"=>181538259, "6"=>170805979, 
					 "7"=>159345973, "8"=>145138636, "9"=>138394717, "10"=>133797422, "11"=>135086622, "12"=>13327530, 
					 "13"=>114364328, "14"=>107043718, "15"=>101991189, "16"=>90338345, "17"=>83257441, "18"=>80373285, 
					 "19"=>58617616, "20"=>64444167, "21"=>46709983, "22"=>50818468, "X"=>156040895, "Y"=>57227415, "MT"=>16569);
} else {
	%chrom_length = ("1"=>249250621, "2"=>243199373, "3"=>198022430, "4"=>191154276, "5"=>180915260, "6"=>171115067, 
					 "7"=>159138663, "8"=>146364022, "9"=>141213431, "10"=>135534747, "11"=>135006516, "12"=>133851895, 
					 "13"=>115169878, "14"=>107349540, "15"=>102531392, "16"=>90354753, "17"=>81195210, "18"=>78077248, 
					 "19"=>59128983, "20"=>63025520, "21"=>48129895, "22"=>51304566, "X"=>155270560, "Y"=>59534049, "MT"=>16569);
}

open (IN, $in_file) or die "open failed for reading: $!";
open (OUT, ">$out_file") or die "open failed for writing: $!";

my $prev_start = 0;
my $prev_stop  = 0;

while (<IN>) {
	chomp;
	my @items = split "\t";
    $prev_stop = $items[1];

	print OUT "$chrom_id\t$prev_start\t$prev_stop\t0\n";
    print OUT "$_\n";

    $prev_start = $items[2];
}

# need to output the last part of the chromosome
#
$prev_stop = $chrom_length{$chrom_id};
print OUT "$chrom_id\t$prev_start\t$prev_stop\t0\n";

close(IN);
close(OUT);

exit;
