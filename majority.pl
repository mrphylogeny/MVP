# majority.pl
# This program demonstrates the use of MVP.
# The arguments are input as options, although they can be
# directly input to subroutines.

#!/usr/bin/perl

use strict;
use warnings;
use Time::HiRes qw ( time );
use Bio::TreeIO;
use Bio::Tree::NodeI;
use Bio::Tree::TreeFunctionsI;
use Getopt::Long;
use MVP;

my $start_run = time();
# default values of options
my $majority_threshold = 51;
my $boot_threshold = 70;
my $feature_selector = 1;
my @features = ();	# feature values of interest

# option definitions
# --maj: feature majority percentage (float)
# --boot: bootstrap value percentage (float)
# --feature: feature code (integer)
# --values: desired values of the feature type (comma-separated values)
my $usage = "Usage: perl $0 [--maj FLOAT1] [--boot FLOAT2] \
			--feature INT --values <comma-separated values> \
			<newick_tree_file>\n";

my $result = GetOptions(
				"maj=f" => \$majority_threshold,
				"boot=f" => \$boot_threshold,
				"feature=i" => \$feature_selector,
				"values=s@" => \&option_handler,
			) or die $usage;

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
}


# main program

# core subroutine
core($tree, $feature_selector);
# utility subroutines
write_collected_node_ids($tree);
write_majority_feature_values($tree, $feature_selector);
write_feature_summary($tree);
write_output_nexus_tree($tree, $boot_threshold, $feature_selector, @features);
find_selected_majority_nodes($tree, $majority_threshold, $feature_selector, @features);
write_tree_summary($tree);

my $end_run = time();
print "Program elapsed time in seconds: ", $end_run - $start_run, "\n";
exit;
