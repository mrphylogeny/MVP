# majority.pl
# This program demonstrates the use of MrP.
# The arguments are input as options, although they can be
# directly input to subroutines.

#!/usr/bin/perl

use strict;
use warnings;
use Bio::TreeIO;
use Bio::Tree::NodeI;
use Bio::Tree::TreeFunctionsI;
use Getopt::Long;
use MrP;

# default values of options
my $majority_threshold = 51;
my $boot_threshold = 70;
my $feature_selector;
my @features = ();	# feature values of interest
my $usage = "Usage: perl $0 [--maj FLOAT1] [--boot FLOAT2] \
			[--date INT1 INT2 | \
			--pri VALUE[,VALUE[,...]] | \
			--sec VALUE[,VALUE[,...]]] \
			<newick_tree_file>\n";

my $result = GetOptions("maj=f" => \$majority_threshold,
				"boot=f" => \$boot_threshold,
				# start year and end year in 4 digits
				"date=i{2}" => \&option_handler,
				"pri=s@" => \&option_handler,
				"sec=s@" => \&option_handler)
				or die $usage;

# check against insufficient command-line arguments
die $usage if (@ARGV < 1);

my $infile = $ARGV[0];
my $intree = Bio::TreeIO->new(-format => "newick", -file => $infile, -internal_node_id => "bootstrap");
# Read the Newick tree in the input file, assuming that there is only one tree.
my $tree = $intree->next_tree;

		
sub option_handler {
	my ($opt_name, @opt_values) = @_;
	
	# Use comma as the option value delimiter.
	# The user can input multiple values concatenated by comma for a
	# multi-value option, e.g. --pri China,USA. The values of a multi-value
	# option are stored in a scalar with a delimiter to separate the elements.
	# To allow an option to take multiple values, a delimiter must be
	# specified. In this case, the comma is used. Therefore the option values
	# are concatenated with comma first, and then split by comma.
	my @values = split /,/, join(",", @opt_values);
	push @features, @values;
	
	# Assign the feature selector with the appropriate value according to
	# the value of the option name.
	if ($opt_name eq "date") {
		$feature_selector = DATE;
	} elsif ($opt_name eq "pri") {
		$feature_selector = PRIMARY;
	} elsif ($opt_name eq "sec") {
		$feature_selector = SECONDARY;
	} else {
		die "Error in option_handler(): unknown option: $!\n";
	}
}


# main program
initialize($tree);
traverse_tree($tree);
write_collected_node_ids($tree);
find_feature_majority($tree, $feature_selector);
write_feature_majority($tree, $feature_selector);
write_feature_summary($tree);
write_output_nexus_tree($tree, $boot_threshold, $feature_selector, @features);
find_selected_majority_nodes($tree, $majority_threshold, $feature_selector, @features);
write_tree_summary($tree);

exit;

