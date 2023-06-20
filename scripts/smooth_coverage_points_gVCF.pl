#!/hgsc_software/perl/perl-5.18.2/bin/perl
#
use lib '/hgsc_software/perl/perl-5.18.2/lib/site_perl/5.18.2/x86_64-linux-thread-multi';
use strict;
use Data::Dumper;

# this script is used to pre-process the binned coverage data for smooth
# the gVCG strategy is used for the calculation 
#

my $in_file  = shift || die "Please specify the input between file\n";
my $out_file = shift || die "Please specify the output file name after smoothing\n";
my $boundary_point = shift || die "Please specify the upper bound used for ExCID run\n";
my $gVCF_PCT   = shift || die "Please specify the percentage used for gVCF formula\n";
my $gVCF_Point = shift || die "Please specify the max point allowed between min and max values\n";

my $upper_bound = $boundary_point;
my $lower_bound = $boundary_point;
my $lowest_allowed = $boundary_point/$gVCF_PCT;
if ($lowest_allowed < 70) { $lowest_allowed = 70; }

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
my $on_flag = 1;
my $prev_chrom_id = 1;

# the field order: 
# chrom_id	start	end		length		coverage
#
open (IN, $in_file) or die "open failed for reading: $!";
open my $out, ">$out_file" or die "open failed for writing: $!";

while (<IN>) {
	chomp;
	next if $_=~/#/;
	my @items = split "\t";

	if ($items[4] > $lowest_allowed) {
		# update the $lower_bound or $upper_bound
		#
		$upper_bound = $items[4] if ($items[4] > $upper_bound);
		$lower_bound = $items[4] if ($items[4] < $lower_bound);

		# check the constraint
		#
		if (($lower_bound*$gVCF_PCT) >= $upper_bound || ($lower_bound + $gVCF_Point) >= $upper_bound) {
			if ($prev_chrom_id eq $items[0]) {
				$on_flag = 2 if $on_flag == 1;

				$blocks{$items[0]}{$items[1]}{$items[2]}{$items[3]} = $items[4];
			} else {
				# a new chrom
				#
				$prev_chrom_id = $items[0];

				# output the previous block info
				#
				if ($on_flag == 2) {
					process_block_info(\%blocks, $out);
				}

				%blocks = ();
				#$lower_bound=$boundary_point;
				#$upper_bound=$boundary_point;

				# now record the new block info
				#
				$blocks{$items[0]}{$items[1]}{$items[2]}{$items[3]} = $items[4];
			}
		} else {
			# doesn't satisfy the BLOCKAVE_Xp constrain
			# break and output block info if block is not null
			#
			if ($on_flag == 2) {
				$on_flag = 1;
				process_block_info(\%blocks, $out);
				%blocks = ();
			}

			$prev_chrom_id = $items[0];
			$lower_bound=$boundary_point;
			$upper_bound=$boundary_point;
			print $out "$_\n";
		}
	} else {
		# fall below the lowest_allowed, so break it
		#
		if ($on_flag == 2) {
			$on_flag = 1;

			process_block_info(\%blocks, $out);
			%blocks = ();
		}

		$prev_chrom_id = $items[0];
		$lower_bound=$boundary_point;
		$upper_bound=$boundary_point;
		print $out "$_\n";
	}
}

close(IN);
close($out);

###################
# sub-routines
###################
sub process_block_info {
	my ($my_block, $output) = @_;
	my ($beginning, $end, $chrom_id);
	my $total_length = 0;
	my $total_coverage = 0;

	foreach my $chr (keys %{$my_block}) {
		$chrom_id = $chr;

		foreach my $start (sort {$a <=> $b} keys %{$my_block->{$chr}}) {
			$beginning = $start if !defined $beginning;

			foreach my $stop (keys %{$my_block->{$chr}{$start}}) {
				$end = $stop;

				foreach my $len (keys %{$my_block->{$chr}{$start}{$stop}}) {
					#$total_length = $total_length + $stop - $start + 1;
					$total_length = $total_length + $stop - $start;
					$total_coverage = $total_coverage + $my_block->{$chr}{$start}{$stop}{$len} * $len;
				}
			}
		}
	}

	if ($total_length == 0) {
		print "stop\n";
	}
	my $average = int ($total_coverage / $total_length + 0.5);

	print $output "$chrom_id\t$beginning\t$end\t$total_length\t$average\n";
}

exit;
