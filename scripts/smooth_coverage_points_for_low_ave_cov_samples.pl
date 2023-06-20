#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;

# this script is used to pre-process the binned coverage data for smooth
# the reason is that for low average coverage samples, the data is too scattered and
# many small length data will fill the graph and make it un-acceptable
#

my $in_file  = shift || die "Please specify the input between file\n";
my $out_file = shift || die "Please specify the output file name after smoothing\n";

# the keys are:
# 1: chrom_id
# 2: start
# 3: end
# 4: length
# 5: coverage
#
my %blocks;

# check to see if we need to group nearby bases
# 1: no, 2: yes
#
my $prev_chrom_id = "ABCDEFG";

# the field order: 
# chrom_id	start	end		length		coverage
# 1       566219  566225  6       101
# 1       2615689 3418597 802908  5
#
open (IN, $in_file) or die "open failed for reading: $!";
open my $out, ">$out_file" or die "open failed for writing: $!";

my ($start, $end, $prev_start, $prev_end) = (0, 0, 0, 0);
my $cur_total_coverage = 0;

while (<IN>) {
	chomp;
	next if $_=~/#/;
	my @items = split "\t";

	# need to check the chromosome id first, if it is different, we need to switch to restart
	#
	if ($prev_chrom_id ne $items[0]) {
		($start, $end, $prev_start, $prev_end) = (0, 0, 0, 0);

		# now output previous chromosome information
		#
		if (length(keys %blocks) > 1) {
			process_block_info(\%blocks, $out);
		}

		$prev_chrom_id = $items[0];
		%blocks = ();
	}

	if ($items[4] >= 10) {
		# need to do the group smoothing
		#
		if ($start == 0) {
			$start = $items[1];
		}

		$end = $items[2];
		$cur_total_coverage = $cur_total_coverage + $items[4]*$items[3];

	} else {
		# just record the block info for both previous block and current block separately
		#

		# for previous block
		#
		if ($start > 0) {
			my $len = $end - $start;
			my $ave_cov = int ($cur_total_coverage / $len + 0.5);
			#$blocks{$items[0]}{$start}{$end}{$len} = $ave_cov;
			print $out "$items[0]\t$start\t$end\t$len\t$ave_cov\n";
			($start, $end, $prev_start, $prev_end) = (0, 0, 0, 0);
			$cur_total_coverage = 0;
		}

		# for current block
		#
		#$blocks{$items[0]}{$items[1]}{$items[2]}{$items[3]} = $items[4];
		print $out "$items[0]\t$items[1]\t$items[2]\t$items[3]\t$items[4]\n";
	}
}

close(IN);
close($out);

###################
# sub-routines
###################
sub process_block_info {
	my ($my_block, $output) = @_;

	foreach my $chr (sort {$a <=> $b} keys %{$my_block}) {

		foreach my $start (sort {$a <=> $b} keys %{$my_block->{$chr}}) {

			foreach my $stop (keys %{$my_block->{$chr}{$start}}) {

				foreach my $len (keys %{$my_block->{$chr}{$start}{$stop}}) {
					print $output "$chr\t$start\t$stop\t$len\t$my_block->{$chr}{$start}{$stop}{$len}\n";
				}
			}
		}
	}
}

exit;
