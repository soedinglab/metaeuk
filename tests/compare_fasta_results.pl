use strict;
use warnings;

=pod

This script is designed to validate the metaeuk fasta results (before and after grouping)

=cut

my ($in_fasta_all_predictions, $in_fasta_grouped_predictions, $in_correct_fasta_all_predictions, $in_correct_fasta_grouped_predictions) = @ARGV;


# read the expected output #
my %as_should_be_metaek_all_predictions;
my %as_should_be_metaeuk_grouped_predictions;
my $as_should_num_preds = parse_metaeuk_fasta_file($in_correct_fasta_all_predictions, \%as_should_be_metaek_all_predictions);
my $as_should_num_grouped_preds = parse_metaeuk_fasta_file($in_correct_fasta_grouped_predictions, \%as_should_be_metaeuk_grouped_predictions);

# read the current output #
my %metaek_all_predictions;
my %metaeuk_grouped_predictions;
my $num_preds = parse_metaeuk_fasta_file($in_fasta_all_predictions, \%metaek_all_predictions);
my $num_grouped_preds = parse_metaeuk_fasta_file($in_fasta_grouped_predictions, \%metaeuk_grouped_predictions);

# compare total number of all predictions:
if ($as_should_num_preds != $num_preds)
{
	die "Failed at all predictions: number of predictions has changed!\n";
}
# compare total number of grouped predictions:
if ($as_should_num_grouped_preds != $num_grouped_preds)
{
	die "Failed at grouped predictions: number of predictions has changed!\n";
}

# compare all results #
compare_as_should_and_current(\%as_should_be_metaek_all_predictions, \%metaek_all_predictions, "all predictions");

# compare grouped results #
compare_as_should_and_current(\%as_should_be_metaeuk_grouped_predictions, \%metaeuk_grouped_predictions, "grouped predictions");

print "ALL OKAY!\n";

sub compare_as_should_and_current
{
	my ($as_should_metaeuk_preds_to_info_ref, $metaeuk_preds_to_info_ref, $type) = @_;
	
	my %as_should_metaeuk_preds_to_info = %{$as_should_metaeuk_preds_to_info_ref};
	my %metaeuk_preds_to_info = %{$metaeuk_preds_to_info_ref};
	
	# compare number of unique predictions #
	if (scalar(keys %as_should_metaeuk_preds_to_info) != scalar(keys %metaeuk_preds_to_info))
	{
		die "Failed at $type: number of unique predictions has changed!\n";
	}
	
	# compare each prediction #
	foreach my $as_should_pred_target (keys %as_should_metaeuk_preds_to_info)
	{
		# check target name #
		if (! exists $metaeuk_preds_to_info{$as_should_pred_target})
		{
			die "Failed at $type: the target '$as_should_pred_target' was not found in result\n";
		}
		
		# go over all contigs (only one in this case but write generally) #
		foreach my $as_should_contig (keys %{$as_should_metaeuk_preds_to_info{$as_should_pred_target}})
		{
			# check contig name #
			if (! exists $metaeuk_preds_to_info{$as_should_pred_target}->{$as_should_contig})
			{
				die "Failed at $type: the contig '$as_should_contig' was not found in result\n";
			}
			
			# go over all strands (only one in this case but write generally) #
			foreach my $as_should_strand (keys %{$as_should_metaeuk_preds_to_info{$as_should_pred_target}->{$as_should_contig}})
			{
				# check strand #
				if (! exists $metaeuk_preds_to_info{$as_should_pred_target}->{$as_should_contig}->{$as_should_strand})
				{
					die "Failed at $type: the strand '$as_should_strand' was not found in result\n";
				}
				
				# if we are here the TCS predicition exists in both - compare content #
				my %as_expected_prediction = %{$as_should_metaeuk_preds_to_info{$as_should_pred_target}->{$as_should_contig}->{$as_should_strand}};
				my %prediction = %{$metaeuk_preds_to_info{$as_should_pred_target}->{$as_should_contig}->{$as_should_strand}};
				
				if ($as_expected_prediction{num_exons} != $prediction{num_exons})
				{
					die "Warning $type: num_exons changed: '$as_should_pred_target' + '$as_should_contig' + '$as_should_strand': from $as_expected_prediction{num_exons} to $prediction{num_exons}\n";
				}
				
				if ($as_expected_prediction{sum_bit_score} != $prediction{sum_bit_score})
				{
					die "Warning $type: bitscore changed: '$as_should_pred_target' + '$as_should_contig' + '$as_should_strand': from $as_expected_prediction{sum_bit_score} to $prediction{sum_bit_score}\n";
				}
				
				if ($as_expected_prediction{evalue} != $prediction{evalue})
				{
					die "Warning $type: evalue changed: '$as_should_pred_target' + '$as_should_contig' + '$as_should_strand': from $as_expected_prediction{evalue} to $prediction{evalue}\n";
				}
				
				if ($as_expected_prediction{low} != $prediction{low})
				{
					die "Warning $type: low coordinate changed: '$as_should_pred_target' + '$as_should_contig' + '$as_should_strand': from $as_expected_prediction{low} to $prediction{low}\n";
				}
				
				if ($as_expected_prediction{high} != $prediction{high})
				{
					die "Warning $type: high coordinate changed: '$as_should_pred_target' + '$as_should_contig' + '$as_should_strand': from $as_expected_prediction{high} to $prediction{high}\n";
				}
				
				# if we are here the overall prediction is the same, now check each exon #
				my $ind = 0;
				foreach my $as_should_exon_ref (@{$as_expected_prediction{exons}})
				{
					my %as_should_curr_exon_info = %{$as_should_exon_ref};
					my %curr_exon_info = %{$prediction{exons}->[$ind]};
					
					foreach my $exon_field_key (keys %as_should_curr_exon_info)
					{
						if ($as_should_curr_exon_info{$exon_field_key} != $curr_exon_info{$exon_field_key})
						{
							die "Warning $type: exon\_$ind changed: '$as_should_pred_target' + '$as_should_contig' + '$as_should_strand': $exon_field_key from $as_should_curr_exon_info{$exon_field_key} to $curr_exon_info{$exon_field_key}\n";
						}
					}
					$ind++;
				}
			}
		}	
	}
}




sub parse_metaeuk_fasta_file
{
	my ($in_metaeuk_united_exons_fasta_file, $metaeuk_preds_to_info_ref) = @_;
	open (my $in, "<", $in_metaeuk_united_exons_fasta_file) or die "could not open $in_metaeuk_united_exons_fasta_file for reading";
	
	my $total_num_preds_in_file = 0;

	my $line = <$in>;
	while (defined $line)
	{
		chomp ($line);
		if ($line =~ m/^>(.*)/)
		{
			$total_num_preds_in_file++;

			# identical_exon1_and_exon4|my_multi_exon_contig|+|398|2.94007e-117|2|100|1444|100[100]:429[429]:330[330]|1184[1202]:1444[1444]:261[243]
			my ($target_id, $contig_id, $strand, $sum_bit_score, $evalue, $num_exons, $low_coord, $high_coord, @exons_info) = split(/\|/, $1);
			my $rounded_evalue = sprintf("%.2g", $evalue);
			$metaeuk_preds_to_info_ref->{$target_id}{$contig_id}{$strand}{sum_bit_score} = $sum_bit_score;
			$metaeuk_preds_to_info_ref->{$target_id}{$contig_id}{$strand}{evalue} = $rounded_evalue;
			$metaeuk_preds_to_info_ref->{$target_id}{$contig_id}{$strand}{num_exons} = $num_exons;
			$metaeuk_preds_to_info_ref->{$target_id}{$contig_id}{$strand}{low} = $low_coord; # MMSeqs2 starts the count at 0, not 1
			$metaeuk_preds_to_info_ref->{$target_id}{$contig_id}{$strand}{high} = $high_coord; # MMSeqs2 starts the count at 0, not 1	
			
			my $exon_num = 0;
			foreach my $exon (@exons_info)
			{
				# 100[100]:429[429]:330[330]
				my ($potential_start, $taken_start, $potential_end, $actual_end, $potential_length, $actual_length);
				if ($exon =~ m/(\d+)\[(\d+)\]\:(\d+)\[(\d+)\]\:(\d+)\[(\d+)\]/)
				{
					$metaeuk_preds_to_info_ref->{$target_id}{$contig_id}{$strand}{exons}[$exon_num] = {potential_start => $1, taken_start => $2, potential_end => $3, actual_end => $4, potential_length => $5, actual_length => $6};
				}
				else
				{
					die "it seems exon_field '$exon' does not match the pattern '\\d\+\[\\d\+\]\:\\d\+\[\\d\+\]\:\\d\+\[\\d\+\]'\n";
				}
				$exon_num++;
			}
		}
		$line = <$in>;
	}
	close ($in);

	return($total_num_preds_in_file);
}