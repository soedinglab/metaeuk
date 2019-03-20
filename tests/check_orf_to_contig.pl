use strict;

use warnings;



my ($in_orfs_to_contigs_tsv, $in_correct_orfs_to_contigs) = @ARGV;

my %orfs_to_contigs;
my %orfs_to_contigs_as_should;

my $raw_num_orfs = parse_orfs_to_contig($in_orfs_to_contigs_tsv, \%orfs_to_contigs);
my $raw_num_orfs_as_should = parse_orfs_to_contig($in_correct_orfs_to_contigs, \%orfs_to_contigs_as_should);



if ($raw_num_orfs != $raw_num_orfs_as_should)
{
	die "not the same number of non-unique records!";
}



if (scalar(keys %orfs_to_contigs) != scalar(keys %orfs_to_contigs_as_should))
{
	die "not the same number of unique records!";
}

compare_hash_to_hash(\%orfs_to_contigs_as_should, \%orfs_to_contigs);

print "ALL OKAY!\n";


sub compare_hash_to_hash
{

	my ($orfs_to_conts_hash1_ref, $orfs_to_conts_hash2_ref) = @_;

	my %orfs_to_conts_hash1 = %{$orfs_to_conts_hash1_ref};
	my %orfs_to_conts_hash2 = %{$orfs_to_conts_hash2_ref};

	foreach my $key1 (keys %orfs_to_conts_hash1)
	{
		if (! exists $orfs_to_conts_hash2{$key1})
		{
			die "$key1 does not exist in target hash!\n";
		}
	}
}



sub parse_orfs_to_contig

{
	my ($in_orfs_to_contigs_file, $orfs_to_conts_hash_ref) = @_;
	
	# can deal with two formats:
	# 0	1	1.000	0.000E+00	0	71	72	2	73	1494	0M
	# 0	2+71	71	my_multi_exon_contig	1	1.000	0.000E+00	0	71	72	2	73	1494	0M

	open (my $in, "<", $in_orfs_to_contigs_file) or die "could not open $in_orfs_to_contigs_file for reading";

	my $line = <$in>;
	my $num_records = 0;
	while (defined $line)
	{
		$line =~ s/\x00//; # remove null byte
		chomp($line);
		if ($line ne '')
		{
			my @line_parts = split(/\s+/, $line);
			
			my $num_elements = @line_parts;
			
			if ($line_parts[$num_elements - 1] eq "0M")
			{
				pop @line_parts;
				$num_elements--;
			}
			
			# only these matter: 0	71	72	2	73	1494
			my @line_parts_that_matter = @line_parts[($num_elements - 6)..($num_elements - 1)];

			$line = join("\t", @line_parts_that_matter);
			$orfs_to_conts_hash_ref->{$line} = 1;
			$num_records++;
		}
		$line = <$in>;
	}

	close ($in);
	return ($num_records);
}

