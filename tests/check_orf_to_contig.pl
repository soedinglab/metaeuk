use strict;
use warnings;

my ($in_dir_to_search_orfs_to_contigs, $in_correct_orfs_to_contigs) = @ARGV;

my ($in_orfs_to_contigs) = glob($in_dir_to_search_orfs_to_contigs . '/*/nucl_6f_orf_aligned_to_contig');

my %orfs_to_contigs;
my %orfs_to_contigs_as_should;

parse_orfs_to_contig($in_orfs_to_contigs, \%orfs_to_contigs);
parse_orfs_to_contig($in_correct_orfs_to_contigs, \%orfs_to_contigs_as_should);

if (scalar(keys %orfs_to_contigs) != scalar(keys %orfs_to_contigs_as_should))
{
	die "not the same number of records!";
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
	
	open (my $in, "<", $in_orfs_to_contigs_file) or die "could not open $in_orfs_to_contigs_file for reading";
	
	my $line = <$in>;
	my $num_records = 0;
	while (defined $line)
	{
		$line =~ s/\x00//; # remove null byte
		chomp($line);
		if ($line ne '')
		{
			$orfs_to_conts_hash_ref->{$line} = 1;
			$num_records++;
		}
		$line = <$in>;
	}
	
	print "\ninserted $num_records in total\n\n";
	close ($in);
}