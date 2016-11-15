# MVP.pm

=pod

=head1 NAME

MVP.pm - Perl package that uses simple majority rule and maximum parsimony
to determine the majority feature values of all nodes of a phylogenetic tree,
given that a feature is selected.

=head1 SYNOPSIS

    #!/usr/bin/perl

    use strict;
    use warnings;
    use Bio::TreeIO;
    use Bio::Tree::NodeI;
    use Bio::Tree::TreeFunctionsI;
    use Getopt::Long;
    use MVP;

    # default values of options
    my $majority_threshold = 51;
    my $boot_threshold = 70;
    my $feature_selector = 1;
    my @features = ();  # feature values of interest
	
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
					"feature=i" =>\$feature_selector,
					"values=s@" =>\option_handler,
				) or die $usage;

    # check against insufficient command-line arguments
    die $usage if (@ARGV < 1);

    my $infile = $ARGV[0];
    my $intree = Bio::TreeIO->new(-format => "newick", \
        -file => $infile, -internal_node_id => "bootstrap");
    # Read the Newick tree in the input file, assuming that there is
    # only one tree.
    my $tree = $intree->next_tree;

            
    sub option_handler {
        my ($opt_name, @opt_values) = @_;
        
        # Use comma as the option value delimiter.
        # The user can input multiple values concatenated by comma
        # for a multi-value option, e.g. --pri China,USA. The values
        # of a multi-value option are stored in a scalar with a
        # delimiter to separate the elements. To allow an option to
        # take multiple values, a delimiter must be specified.
        # In this case, the comma is used. Therefore the option
        # values are concatenated with comma first, and then split
        # by comma.
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
    write_output_nexus_tree($tree, $boot_threshold, \
        $feature_selector, @features);
    find_selected_majority_nodes($tree, $majority_threshold, \
        $feature_selector, @features);
    write_tree_summary($tree);

    exit;



=head1 DESCRIPTION

This package determines and highlights the values of a selected feature of
a phylogenetic tree. The features are inserted into the node identifiers in
a specific format. The selected values are determined by simple majority
rule. If there is tie in the vote of the feature values, the majority feature
value is determined by maximum parsimony.

A node identifier of a tree can contain values of some characteristics, for
example, geographical region where a taxon is found. These characteristics
are called features, and different kinds of features are inserted into the
node identifiers in a specific format. For a selected feature, we are
interested to know what the majority of a particular node is. This majority
feature value can be determined by counting the feature values of all its
descendant nodes, and taking the value with the highest vote. This majority
determination procedure goes from the leaf nodes to the root.

Sometimes there may be tie in the number of votes at a particular node. For
example, for a feature X, a node N has feature values of one A and one B. To
resolve the tie, the majority feature value of N's ancestor is consulted. If
the value is not available, this process continues recursively until a
majority feature value is obtained. If the majority feature value obtained
is A, N takes this value and hence its majority feature value is now A. Since
this procedure results in the fewest changes in majority feature values, tie
votes are resolved by maximum parsimony. If, however, the ancestor's majority
feature value is neither A nor B, or undetermined, the ancestor of N's ancestor
is consulted for its majority feature value, until a value is obtained.

In a rare case, a node may refer to its ancestors recursively up to the root
and still have the majority feature value undetermined. If this happens, the
value "NA" is assigned as the majority feature value of the root, and
subsequently down to all nodes whose majority feature value are previously
undetermined.

The majority determination procedure consists of three steps, which
corresponds to the three core subroutines as described in the SUBROUTINES
section.

=over 3

=item Initialization

Each node is assigned a unique node number, and variables are initialized.
The assignment of node numbers is in a pre-order manner.

=item Tree traversal

The tree is traversed from the root to the leaf nodes in a depth-first,
pre-order manner. For each node, the node identifiers of its descendant nodes
are collected. This process goes from the leaves back to the root.

=item Majority feature value determination

For each node, the values of a selected feature are extracted and counted.
The value with the highest count, by simple majority rule, is regarded as
the majority feature value. If tie votes are present, the node's ancestor
is consulted until a majority feature value is obtained, and this value is
assigned as the feature majority value of that node.

=back

After majority feature value determination, users can obtain text reports
and graphic visualization with the optional utility subroutines.

=head1 INPUT

The input requires a phylogenetic tree file in Newick format. Here the leaf
node identifiers are composed of fields. Each field represents a feature
which contains a value of that particular node. The node identifier is in
the format

accession_number@date@feature1@feature2@feature3@...

where the feature values are separated by the field delimiter @,
for example, AB1234@1997/06/30@Hong_Kong@Human.

The accession number and date are mandatory, which can follow with an
arbitrary number of optional features. The date is in yyyy/mm/dd or,
at least, yyyy format, in which yyyy (year) has 4 digits, mm (month) 2 digits,
and dd (day) 2 digits. The @ sign is used because it seldoms appears in the
content of a Newick tree, and it is a valid character in Newick format.
It should be noted that a field delimiter cannot be a word separator in a
node identifier simultaneously.
<<<<<<< 84498933f8834c1e694a2fb8f0201f443907aa33

Each field of a node identifier can be selected with an integer, which
starts with zero. For a node identifier with N fields, the first feature
is 0, and the last feature is N - 1. Therefore 0 represents accession number,
and 1 represents date, and so on.
=======
>>>>>>> Renamed write_feature_majority to write_majority_feature_values()

Each field of a node identifier can be selected with an integer, which
starts with zero. For a node identifier with N fields, the first feature
is 0, and the last feature is N - 1. Therefore 0 represents accession number,
and 1 represents date, and so on.

=head1 SUBROUTINES

The package contains three core subroutines to be run sequentially as follows.
One feature can be selected at a time.

=head2 CORE SUBROUTINES

=over

=item core($tree, $feature_selector)

$tree: input Newick tree

$feature_selector: integer, default 1 (date)

Functions: It includes the core functions of the package: initialization,
tree traversal, and majority feature value determination.

=back

The core() subroutine runs these three subroutines sequentially.

=over

=item _initialize($tree)

$tree: input Newick tree

Functions: It assigns a unique node number to each node and initializes
variables and tag values of each node.

=item _traverse_tree($tree)

$tree: input Newick tree

Functions: It traverses the tree from the root to the leaves, and collect
node identifiers from the leaves to the root. At each node, its descendant
nodes' identifiers are concatenated by comma to form a string.

=item _find_majority_feature_values($tree, $feature_selector)

$tree: input Newick tree

$feature_selector: integer, default 1 (date)

Functions: It extracts the values of the selected feature and determines the
majority feature value of each node by simple majority rule and maximum
parsimony. The feature selector corresponds to the feature fields in a node
identifier, starting with zero. Therefore 0 represents accession number, and
1 represents date, and so on.

=back

=head2 UTILITY SUBROUTINES

Utility subroutines are optional subroutines which generate text reports and
graphical visualization for users. Since they depend on the majority feature
values to produce output, they must be run after running the three core
subroutines.

=over

=item write_output_nexus_tree($tree, $boot_threshold, $feature_selector,
@targets)

$tree: input Newick tree

$boot_threshold: bootstrap value threshold, default 70

$feature_selector: integer, default 1 (date)

@targets: target values of a selected feature (default ALL), which are
supplied as a list of comma-separated values

Functions: It writes an output tree in NEXUS format which is ready for
visualization with FigTree (http://tree.bio.ed.ac.uk/software/figtree/).
Selected nodes are highlighted. If the majority feature values of the nodes
match the user's target values, the branches are coloured in green with the
majority feature percentages shown on the branches. If empty target values
are given, it is assumed that all majority feature values are selected. Among
those selected branches, if the nodes have bootstrap values greater than or
equal to the threshold, the branches and the bootstrap values are coloured
in red. Unselected branches are shown in black with their majority feature
values masked as "NA".

Output: outtree.nexus

=item write_majority_feature_values($tree, $feature_selector)

$tree: input Newick tree

$feature_selector: integer, default 1 (date)

Functions: It writes the majority feature values and their percentages of the
selected feature. If the majority feature value is also the absolute majority
feature value (>50% of the collected feature values of a node), it is also
shown in the absolute majority column; otherwise "nil" is shown. In addition,
the immediate ancestor of each node is shown.

Output: majority_feature.txt

=item find_selected_majority_nodes($tree, $majority_threshold,
$feature_selector, @targets)

$tree: input Newick tree

$majority_threshold: majority feature percentage threshold, default 51

$feature_selector: integer, default 1 (date)

@targets: target values of a selected feature (default ALL), which are 
supplied as a list of comma-separated values

Functions: It writes the nodes whose majority feature values of the selected
feature match the user's target values and have the majority feature percentage
greater than or equal to the threshold. If empty target values are given, it
is assumed that all majority feature values are selected. For each selected
node, its descendant nodes are shown. If a selected node has no descendants,
"nil" is shown.

Output: selected_nodes.txt

=item write_collected_node_ids($tree)

$tree: input Newick tree

Functions: It writes the collected node identifiers of all nodes after
traversal. Collected node identifiers are concatenated by comma. The
identifiers of internal nodes are removed after traversal.

Output: collected_node_ids.txt

=item write_feature_summary($tree)

$tree: input Newick tree

Functions: It writes all values and their counts for each feature in the tree,
and sums a subtotal count of that feature.

Output: feature_summary.txt

=item get_non_binary_nodes($tree)

$tree: input Newick tree

Functions: It reads the input tree and finds whether a node is binary, i.e.
bifurcated to two branches. If a node has more than two branches, it is stored
in a list and the whole list is returned.

Return type: a list of Node objects

=item write_tree_summary($tree)

$tree: input Newick tree

Functions: It prints a brief summary of the tree on the screen, including
total number of nodes, number of leaf nodes, tree height, and total branch
length. If the tree is not binary (at least one node has three branches or
more), the non-binary nodes and their descendant nodes are also shown.

=back

=head1 LOGS

During the majority feature value determination procedure, intermediate data
are recorded in logs, which are mainly for developers' use for information
and diagnostics.

=over

=item traversal_log.txt

It records how the input tree is traversed. It shows how the nodes are
visited and how relevant tag values of nodes are changed.

=item collected_node_ids.txt

It records what node identifiers are collected at each node after travsersal.
Collected node identifiers are concatenated by comma. The identifiers of
internal nodes are removed after traversal.

=item feature_counts.txt

It lists the raw data produced in the majority feature value determination
procedure at each node. It also records how ties in feature values are
resolved and which ancestors are referred to.

=back

=head1 DEPENDENCIES

Math::Round is required. GetOpt::Long is recommended to pass arguments
to subroutines although it is not mandatory.

=head1 AUTHOR

Chi Ho WONG (genewch@yahoo.com)

=cut


#!/usr/bin/perl

package MVP;

use Math::Round;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
			core
            get_non_binary_nodes
			write_collected_node_ids
            find_feature_majority
            write_majority_feature_values
            write_feature_summary
            find_selected_majority_nodes
            write_output_nexus_tree
            write_tree_summary
			);

use strict;

# Global variable that stores the maximum number of feature fields found
our $max_fields = 0;

# Run the core subroutines sequentially.
sub core {
	my ($tree, $feature_selector) = @_;
	_initialize($tree);
	$max_fields = _get_max_fields($tree);
	_traverse_tree($tree);
	_find_majority_feature_values($tree, $feature_selector);
}


# Initialize tags of an input tree.
sub _initialize {
    my ($tree) = @_;
    my @nodes = $tree->get_nodes;
    my $root = $tree->get_root_node;
    my $i = 0;  # node number, zero as the root
    
    print "Initializing parameters for nodes and assigning node numbers.\n";
    # assign node number, and initialize tags
    # traversal is depth-first, pre-order traversal by default in BioPerl
    foreach my $node (@nodes) {
        # initialize tags with default values
        $node->add_tag_value("node_number", $i);
        $node->add_tag_value("as_leaf", 0);
        $node->add_tag_value("ties_in_hi_count", 0);
        $node->add_tag_value("majority", "none");
        $node->add_tag_value("majority_percentage", 0);
        $node->add_tag_value("absolute_majority", "nil");
        $node->add_tag_value("hi_count", 0);
        $node->add_tag_value("total_count", 0);
        # determine node type which is leaf, internal, or root
        if ($node->is_Leaf) {
            $node->add_tag_value("node_type", "leaf");
            $node->set_tag_value("as_leaf", 1); # true for a leaf
            # set node id as the default tag value of a leaf node
            $node->add_tag_value("collected_node_ids", $node->id);
        } elsif ($node eq $root) {
            $node->id($i);
            $node->add_tag_value("node_type", "root");
            # define the initial node identifier as node_number@init
            $node->add_tag_value("collected_node_ids", $node->get_tag_values("node_number") . "\@init");
        }else {
            $node->id($i);
            $node->add_tag_value("node_type", "internal");
            # define the initial node identifier as node_number@init
            $node->add_tag_value("collected_node_ids", $node->get_tag_values("node_number") . "\@init");
        }
        $i++;
    }
}


# Return a list of nodes which are not bifurcated.
sub get_non_binary_nodes {
    my ($tree) = @_;
    my @non_binary_nodes = ();  # list of non-binary nodes
    my $all_binary = 1; # flag whether all nodes are binary
    my @nodes = $tree->get_nodes;
    
    foreach my $node (@nodes) {
        # check whether a node in a tree is binary
        if (!$tree->is_binary($node)) {
            $all_binary = 0;
            push @non_binary_nodes, $node;
        }
    }
    return @non_binary_nodes;
}


# Return the maximum number of feature fields found.
sub _get_max_fields {
	my ($tree) = @_;
	my @leaves = $tree->get_leaf_nodes;
	my $max = 0;	# maximum number of fields in a leaf id
	
	foreach my $leaf (@leaves) {
		# If a leaf id does not contain the field delimiter,
		# print an error message and terminate.
		if ($leaf->id !~ /@/) {
			die "Error in _get_max_fields(): field delimiter in leaf node identifiers not found.\n";
		} else {
			my $n = scalar split /@/, $leaf->id;
			if ($n > $max) {
				$max = $n;
			}
		}
	}
	return $max;
}


# The subroutine that controls _traverse().
sub _traverse_tree {
    my ($tree) = @_;
    my $root = $tree->get_root_node;
    my $traversal_log = "traversal_log.txt";    # traversal log
	
    print "Writing traversal log to $traversal_log.\n";
    open my $fh, ">", $traversal_log or die $!;
	print $fh "#Traversal log: $traversal_log\n";
	print $fh "#Input file: $ARGV[0]\n";
	print $fh "#", "=" x 40, "\n";
    # Traverse the tree from the root to the leaves.
    # Traversal is depth-first, pre-order by default in BioPerl.
    _traverse($fh, $root);
    close $fh;
}


# Traverse the input tree.
sub _traverse {
    my ($fh, $start_node) = @_;
	my @descendants = $start_node->each_Descendent;
	# Flag whether all descendants of a start node are leaves
    my $all_descendents_as_leaves = 1;
    
    # diagnostics: print the start node and whether it can be
	# treated as a leaf node
    print $fh "Start node in _traverse() is ", $start_node->get_tag_values("node_number"), ", as_leaf = ", $start_node->get_tag_values("as_leaf"), "\n";
    
    foreach my $node (@descendants) {
        print $fh "\t", "Go to node ", $node->get_tag_values("node_number"), ", as_leaf = ", $node->get_tag_values("as_leaf"), "\n";
        
		# If $node is not a leaf node, traverse recursively.
        if ($node->get_tag_values("as_leaf") == 0) {
            _traverse($fh, $node);
        }
        print $fh "\t\t", "===>Exit from node ", $node->get_tag_values("node_number"), ", as_leaf = ", $node->get_tag_values("as_leaf"), ", collected node ids: ", $node->get_tag_values("collected_node_ids"), "\n";
    }
    # If all descendants are treated as leaves, set the start node
	# as a leaf so that the tree can be traversed back to the root.
    ALL_ARE_LEAVES: foreach my $node (@descendants) {
        if ($node->get_tag_values("as_leaf") == 0) {
            $all_descendents_as_leaves = 0;
            last ALL_ARE_LEAVES;
        }
        # treat the start node as a leaf
        $start_node->set_tag_value("as_leaf", 1);
    }
    # If all descendants of $node are regarded as leaves,
	# add all descendants' node identifiers to $node.
    if ($all_descendents_as_leaves == 1) {
        foreach my $descendant (@descendants) {
            _add_node_id($fh, $start_node, $descendant);
        }
    }
}


# Add the node identifiers of a node's descendants to that node.
sub _add_node_id {
    my ($fh, $current_node, $child) = @_;   
    my $current_node_collected_node_ids = $current_node->get_tag_values("collected_node_ids");
    my $child_collected_node_ids = $child->get_tag_values("collected_node_ids");
    
    # Join current node's collected node ids and child's collected node ids
    # with a comma to form a comma-separated value. A delimiter is needed to
    # distinguish among the collected node ids from different nodes.
    my $added_ids = join(",", $current_node_collected_node_ids, $child_collected_node_ids);
    
    # Remove all the default collected_node_ids tags after concatenation to
    # prevent the warning for empty strings from being triggered by
    # get_tag_values().
    # The format of the default collected_node_ids tag value is defined in
    # initialize(). The default tag value cannot be set as all digits only,
    # because an all-digit tag value may cause all digits in the leaf node id
    # to be removed as well. Therefore a non-digit suffix is added to the
    # default tag value.
	
	# The initial identifiers of the internal nodes must be matched first,
	# and then be removed.
    $added_ids =~ s/\d+\@init//g;
    
    # set_tag_value() is needed because a child's collected_node_ids tag
    # cannot be directly added to its ancestor's features tag with
    # add_tag_value(). Therefore the statement
    # $current_node->add_tag_value("collected_node_ids", \
    # $child->get_tag_values("collected_node_ids")); does not work.
    $current_node->set_tag_value("collected_node_ids", $added_ids);
    print $fh "\t", "Node ", $child->get_tag_values("node_number"), "'s ancestor is ", $current_node->get_tag_values("node_number"), " which has node ids: ", $current_node->get_tag_values("collected_node_ids"), "\n";
    
    # store the ancestor of the current node to a temporary variable
    my $temp = $current_node->ancestor;
    # update the current node as a child node
    $child = $current_node;
    # update the current node's ancestor as the current node
    $current_node = $temp;
    # diagnostics: print current node and child after update
    if ($current_node) {
        print $fh "\t", "----current_node is ", $current_node->get_tag_values("node_number"), " child is ", $child->get_tag_values("node_number"), "\n";
    }
}


# Write the collected node identifiers of each node.
sub write_collected_node_ids {
    my ($tree) = @_;
    my @nodes = $tree->get_nodes;
    # file containing collected node ids of all nodes after traversal
    my $collected_node_ids = "collected_node_ids.txt";
    
    print "Writing collected node ids of all nodes to $collected_node_ids.\n";
    open my $fh, ">", $collected_node_ids or die $!;
    # write report header
	print $fh "#Collected node identifiers after traversal: $collected_node_ids\n";
	print $fh "#Input file: $ARGV[0]\n";
	print $fh "#", "=" x 40, "\n";
    print $fh "#node_number", "\t", "node_id", "\t", "type", "\t", "collected_node_ids", "\t", "node_count", "\n";
    foreach my $node (@nodes) {
        print $fh $node->get_tag_values("node_number"), "\t", $node->id, "\t", $node->get_tag_values("node_type"), "\t", $node->get_tag_values("collected_node_ids"), "\t", my $node_count = split (",", $node->get_tag_values("collected_node_ids")), "\n";
    }
    close $fh;
}


# Determine the majority feature value of a selected feature.
sub _find_majority_feature_values {
    my ($tree, $feature_selector) = @_;
    my @nodes = $tree->get_nodes;
    my $root = $tree->get_root_node;
    # file that records the majority feature determination procedure
    my $feature_counts = "feature_counts.txt";
    # default selected feature type is date
    $feature_selector = 1 if ($feature_selector eq "");
    
	# Prevent an invalid feature from being selected.
	if ($feature_selector >= $max_fields || $feature_selector < 0) {
		die "Invalid feature code ($feature_selector) is entered. The feature code exceeds the maximum number of fields found ($max_fields) or is negative.\n";
	}
    
    print "Writing raw feature counts of all nodes to $feature_counts.\n";
    open my $fh, ">", $feature_counts or die $!;
    # write report header
	print $fh "#Feature value counting and majority determination: $feature_counts\n";
	print $fh "#Input file: $ARGV[0]\n";
    print $fh "#Selected feature code: $feature_selector\n";
    print $fh "#", "=" x 40, "\n";
    print $fh "#node_number", "\t", "feature_value", "\n";

    foreach my $node (@nodes) {
        # hash to count each feature with its feature value as the key
        my %counts = ();
        my $sum = 0;    # sum of feature counts
        my @feature_values = ();    # feature values obtained from the node
        
        # get the feature values of a node according to the selected feature
        push @feature_values, _get_feature_values($node, $feature_selector);
        
		# count the number of each feature value
        foreach my $value (@feature_values) {
            # print the node number and its extracted feature value
            print $fh $node->get_tag_values("node_number"), "\t", $value, "\n";
            # increment the counter for the feature name by 1
            $counts{$value}++;
            # add to the sum of feature counts
            $sum++;
        }
        
        my $max = 0;    # maximum feature count
        my @majority = ();  # store the features with the highest count
        
        # sort features by count in descending order, and find the features
        # having the highest counts, i.e. feature majority values
        foreach my $c (sort { $counts{$b} <=> $counts{$a} } keys %counts) {
            # find the highest feature count
            $max = $counts{$c} > $max ? $counts{$c} : $max;
            # store the feature names with the highest counts
            push @majority, $c if ($counts{$c} == $max);
            print $fh $c, " has count ", $counts{$c}, "\n";
            # find the feature value with absolute majority whose count is
            # > 50% of the total count
            $node->set_tag_value("absolute_majority", $c) if ($counts{$c} / $sum > 0.5);
            
            # Set the node's ties_in_hi_count tag to 1 if there are ties in
            # highest feature counts, but not ties in non-highest feature
            # counts, i.e. there are 2 or more feature values having the same
            # highest counts. Therefore @majority should have 2 or more
            # elements.
            if (scalar @majority >= 2) {
                $node->set_tag_value("ties_in_hi_count", 1);
            }
        }
        print $fh "Sum of feature counts = $sum\n";
        
        # assign the feature majority name to the majority tag
        $node->set_tag_value("majority", @majority);
        $node->set_tag_value("hi_count", $max);
        $node->set_tag_value("total_count", $sum);
        
        # assign the node's majority feature percentage to its majority
        # percentage tag
        $node->set_tag_value("majority_percentage", $max / $sum * 100);
        
        # Since a feature with absolute majority has count > 50%, it must not
        # have ties in feature count with other features. Therefore the number
        # of elements of @majority is 1, i.e. there should be only one feature
        # having absolute majority.
        if ($max / $sum > 0.5) {
            print $fh $node->get_tag_values("absolute_majority"), " has absolute majority (>50%) of ", round($node->get_tag_values("majority_percentage")), "% (", $max, "/", $sum, ")\n";
        }
            
        # show that there are ties of highest feature counts, not ties of
        # non-highest feature counts
        if ($node->get_tag_values("ties_in_hi_count") == 1) {
            print $fh "\t", "****There are ties for features with the highest count $max\n";
        }
            
        # print the majority feature value
        print $fh "\t", "####The majority feature value is ";
        for (my $i = 0; $i < scalar @majority; $i++) {
            print $fh $majority[$i];
            # print the comma if there are 2 or more majority feature values
            # and the current feature is not the last element
            print $fh ", " if (scalar @majority >= 2 && $i != $#majority);
        }
        print $fh "\n";
        
        # if there are 2 or more feature majority values, look up the
        # feature majority value of the node's ancestor recursively
        if (scalar @majority >= 2) {
            my $vote = _get_simple_majority_vote($fh, $node, $root);
            # assign the simple majority vote to the majority feature value
            # of the node
            $node->set_tag_value("majority", $vote);
            print $fh "\t", "The simple majority vote result is ", $node->get_tag_values("majority"), "\n";
        }
        
        # write internal node id as node_number@feature_majority
        if (!$node->is_Leaf) {
            $node->id($node->get_tag_values("node_number") . "@" . $node->get_tag_values("majority"));
        }
    }
    close $fh;
}


# Return the majority feature value by simple majority vote.
sub _get_simple_majority_vote {
    my ($fh, $node, $root) = @_;
    
    # This is the base case of recursion.
	# If the node reaches the root and still does not get a majority
	# feature value, return "NA" as the majority feature value.
    if ($node eq $root) {
        return "NA";
    }
	
    # diagnostics: print the node's ancestor and the majority feature value
	# of the ancestor
    print $fh "\t", "Node ", $node->get_tag_values("node_number"), "'s ancestor is ", $node->ancestor->get_tag_values("node_number"), " which has the majority feature value ";
    my @majority = $node->ancestor->get_tag_values("majority");
    for (my $i = 0; $i < scalar @majority; $i++) {
        print $fh $majority[$i];
        # print the comma if there are 2 or more majority feature values and
        # the current feature value is not the last element
        print $fh ", " if (scalar @majority >= 2 && $i != $#majority);
    }
    print $fh "\n";
    
    # if $node's ancestor have two or more feature majority, $node's ancestor
    # goes into _get_simple_majority_vote() recursively until one majority
    # feature is obtained.
    if (scalar (my @m = $node->ancestor->get_tag_values("majority")) >= 2) {
        return _get_simple_majority_vote($fh, $node->ancestor);
    }
    
    # if the node reaches the root and the node has finally no majority
    # feature value, return NA as the majority feature value
    return $node->ancestor->get_tag_values("majority");
}


# Given a selected feature, return a list of feature values of a node.
sub _get_feature_values {
    my ($node, $feature_selector) = @_;
	# The collected node ids are in comma-separated values.
	# Split the collected node ids by comma, and store the split values
	# in @collected_node_ids.
    my @collected_node_ids = split /,/, $node->get_tag_values("collected_node_ids");
	# array that stores the extracted feature values to be returned
    my @returned_features = ();
	
	# The leaf node id is in the format
    # accession_number@yyyy/mm/dd@feature1@feature2@feature3@...
    # The field delimiter is @. The accession number and date are mandatory.
	# Other fields are optional.
    
    foreach my $id (@collected_node_ids) {
        # only process non-blank node identifiers after splitting by comma
        if ($id) {
            # Ensure that the field delimiters have to be found
			if ($id !~ /@/) {
				die "Error in _get_feature_values(): field delimiter not found in node identifiers.\n";
			} else {
				my @split_node_id = split /@/, $id;
				# If a node identifier has fewer fields than the
				# maximum number of fields found, add dummy fields with
				# the value "dummy" up to the maximum number of fields.
				if (scalar @split_node_id < $max_fields) {
					for (my $i = scalar @split_node_id; $i < $max_fields; $i++) {
						$split_node_id[$i] = "dummy";
					}
				}
				# $feature_selector is the subscript of @split_node_id
				# which points to the selected feature.
				# Therefore $split_node_id[0] is accession number,
				# $split_node_id[1] is date,
				# $split_node_id[2] is feature 1,
				# $split_node_id[3] is feature 2, etc.
				
				# If date is selected, the year part is extracted.
				# Otherwise the selected feature is extracted.
				if ($feature_selector == 1) {
					# split the date by the date separator, i.e. /
					my @split_date = split /\//, $split_node_id[1];
					# return the year part of the date, i.e. $split_date[0]
					push @returned_features, $split_date[0];
				} else {
					push @returned_features, $split_node_id[$feature_selector];
				}
				# Therefore $split_node_id[0] is accession number,
				# $split_node_id[1] is date,
				# $split_node_id[2] is feature 1,
				# $split_node_id[3] is feature 2, etc.
				
				# If date is selected, the year part is extracted.
				# Otherwise the selected feature is extracted.
				if ($feature_selector == 1) {
					# split the date by the date separator, i.e. /
					my @split_date = split /\//, $split_node_id[1];
					# return the year part of the date, i.e. $split_date[0]
					push @returned_features, $split_date[0];
				} else {
					push @returned_features, $split_node_id[$feature_selector];
				}
			}
        }
    }
	
    return @returned_features;
}


# Write the majority feature values of a selected feature.
sub write_majority_feature_values {
    my ($tree, $feature_selector) = @_;
    my @nodes = $tree->get_nodes;
    my $root = $tree->get_root_node;
    my $majority_feature = "majority_feature.txt";
    # default selected feature type is date
    $feature_selector = 1 if ($feature_selector eq "");
    
	print "Selected feature code: $feature_selector\n";
    print "Writing majority feature values to $majority_feature.\n";
    open my $fh, ">", $majority_feature or die $!;
    # write report header
	print $fh "#Majority feature value report: $majority_feature\n";
	print $fh "#Input file: $ARGV[0]\n";
    print $fh "#Selected feature code: $feature_selector\n";
    print $fh "#", "=" x 40, "\n";
    print $fh "#node_number", "\t", "node_id", "\t", "ancestor", "\t", "majority", "\t", "count", "\t", "total", "\t", "majority%", "\t", "absolute_majority_(>50%)", "\n";
    
    foreach my $node (@nodes) {
        # Since the root has no ancestor, its ancestor is printed as "no_ancestor".
        print $fh $node->get_tag_values("node_number"), "\t", $node->id, "\t", ($node != $root) ? $node->ancestor->id : "no_ancestor", "\t", $node->get_tag_values("majority"), "\t", $node->get_tag_values("hi_count"), "\t", $node->get_tag_values("total_count"), "\t", round($node->get_tag_values("majority_percentage")), "\t", $node->get_tag_values("absolute_majority"), "\n";
    }
    close $fh;
}


# Write a summary for all features.
sub write_feature_summary {
    my ($tree) = @_;
    my @leaves = $tree->get_leaf_nodes;
    my $feature_summary = "feature_summary.txt";
    
    print "Writing feature summary to $feature_summary.\n";
    open my $fh, ">", $feature_summary or die $!;
	# Write report header
	print $fh "#Feature summary: $feature_summary\n";
	print $fh "#Input file: $ARGV[0]\n";
	print $fh "#Maximum number of feature fields: $max_fields\n";
    print $fh "#", "=" x 40, "\n";
	
    for (my $code = 0; $code < $max_fields; $code++) {
        my @features = ();
		# hash that counts the number of different feature values
        my %all_features = ();
        
        # get the feature values of each leaf node
        foreach my $leaf (@leaves) {
            push @features, _get_feature_values($leaf, $code);
        }
        
        # count the frequency of each feature value
        foreach my $feature (@features) {
            $all_features{$feature}++;
        }
        
        # count the number of different feature values
        my $count = 0;
        foreach my $key (keys %all_features) {
            $count++;
        }
        
        print $fh "#Feature code: $code\n";
        print $fh "#Number of features: $count\n";
        print $fh "#feature_value", "\t", "count", "\n";
        foreach my $key (sort {$a cmp $b} keys %all_features) {
            print $fh $key, "\t", $all_features{$key}, "\n";
        }
        print $fh "#", "-" x 40, "\n";
    }
    close $fh;
}


# Given a selected feature, a majority feature percentage, and target values,
# find the nodes whose majority feature percentage is greater than or
# equal to the threshold and whose majority feature value matches the
# target values.
sub find_selected_majority_nodes {
    my ($tree, $majority_percentage, $feature_selector, @targets) = @_;
    my @nodes = $tree->get_nodes;
    my @selected_nodes = ();
    my @final_nodes = ();
    my $selected_nodes = "selected_nodes.txt";
    # default majority percentage is 51
    $majority_percentage = 51 if (!$majority_percentage);
    # default selected feature type is date
    $feature_selector = 1 if ($feature_selector eq "");
    # default target feature value is ALL, i.e. all values are selected
    @targets = ("ALL") if (! @targets);
    
	# Check arguments if date is selected to determine feature majority.
	_check_date_arguments($feature_selector, @targets);
	
	# Select nodes whose majority percentage is greater than or equal to the threshold.
    @selected_nodes = grep {$_->get_tag_values("majority_percentage") >= $majority_percentage} @nodes;
    
    # If the selected feature type is not date, further select the nodes
    # in which the majority feature value matches the target feature values.
    # If the selected feature type is date, further select the nodes in which
    # the majority feature is within the date range.
    foreach my $s (@selected_nodes) {
        if ($feature_selector != 1) {
            MATCH_SELECTED_NODES: foreach my $target (@targets) {
                if ($target eq "ALL" || $s->get_tag_values("majority") =~ /$target/i) {
                    push @final_nodes, $s;
                    last MATCH_SELECTED_NODES;
                }
            }
        } else {
			# For dates, only the first two values of @targets are used.
			if ($s->get_tag_values("majority") >= $targets[0] && $s->get_tag_values("majority") <= $targets[1]) {
                    push @final_nodes, $s;
            }
        }
    }
    
    print "Writing selected nodes with majority feature value(s) = @targets and majority % >= $majority_percentage to $selected_nodes.\n";
    open my $fh, ">", $selected_nodes or die $!;
    # write report header
	print $fh "#Selected feature values report: $selected_nodes\n";
	print $fh "#Input file: $ARGV[0]\n";
    print $fh "#Threshold majority %: $majority_percentage\n";
    print $fh "#Selected feature code: $feature_selector\n";
    print $fh "#Selected feature value(s):  @targets\n";
    print $fh "#Number of selected nodes: ", scalar @final_nodes, "\n";
	print $fh "#", "=" x 40, "\n";
    print $fh "#node_number", "\t", "node_id", "\t", "majority", "\t", "majority%", "\t", "descendants", "\n";
    
    # write all descendants of the selected nodes in a "flattened" format
    # in parentheses
    foreach my $unit (@final_nodes) {
        print $fh $unit->get_tag_values("node_number"), "\t", $unit->id, "\t", $unit->get_tag_values("majority"), "\t", round($unit->get_tag_values("majority_percentage")), "\t";
        print $fh "(";
        # flag whether the first node id is printed, initialized as false
        my $first = 0;
        print $fh "nil" if ($unit->is_Leaf);
        foreach my $subnode ($unit->get_all_Descendents()) {
            print $fh "," if ($first == 1);
            $first = 1; # set as true when the first node id is printed              
            print $fh $subnode->id;
        }
        print $fh ")\n";
    }
}


# Given a selected feature and a bootstrap value threshold,
# write an output tree in NEXUS format.
sub write_output_nexus_tree {
    my ($tree, $boot_threshold, $feature_selector, @targets) = @_;
    # output NEXUS tree file with parameters for FigTree
    my $outfile = "outtree.nexus";
    
    # Make a copy of the original tree and work on this copy to prevent
	# original data from being overwritten.
    my $tree_copy = $tree->clone;
    my @nodes = $tree_copy->get_nodes;
    my @leaves = $tree_copy->get_leaf_nodes;
    # default bootstrap value threshold is 70
    $boot_threshold = 70 if (!$boot_threshold);
    # default selected feature type is date
    $feature_selector = 1 if ($feature_selector eq "");
    # default target feature value is ALL, i.e. all values are selected
    @targets = ("ALL") if (! @targets);
    
	# Check arguments if date is selected to determine feature majority.
	_check_date_arguments($feature_selector, @targets);
	
    print "Writing NEXUS output tree with bootstrap value >= $boot_threshold and selected feature value(s) = @targets to file $outfile.\n";
    open my $fh, ">", $outfile or die $!;
    # Part 1: write the TAXA block
    print $fh "#NEXUS\n";
    print $fh "BEGIN TAXA;\n";
    print $fh "\t", "DIMENSIONS NTAX=", scalar $tree_copy->get_leaf_nodes, ";\n";
    print $fh "\t", "TAXLABELS\n";
    
    # add colour code after the leaf node id where appropriate
    foreach my $leaf (@leaves) {
        my $colour_attribute = "";
        my $green = "#00ff00";  # green
        
        # If the feature type is not date, prepare the colour attribute for
        # a node with the matching target feature value.
        # If the feature type is date, extract the year from the leaf node id,
        # and then prepare the colour attribute for that node.
        if ($feature_selector != 1) {
            TAXA_MATCH_TARGETS: foreach my $target (@targets) {
				my @split_leaf_id = split /@/, $leaf->id;
                # colour a target feature's leaf node; otherwise no colour
                if ($split_leaf_id[$feature_selector] =~ /$target/i) {
                    $colour_attribute = "[&!color=" . $green . "]"; # green
                    last TAXA_MATCH_TARGETS;
                } else {
                    $colour_attribute = "";
                }
            }
        } else {
            my @split_node_id = split /\@/, $leaf->id;
            # the date is $split_node_id[1], and the year part is $split_date[0]
            my @split_date = split /\//, $split_node_id[1];
			# For date, only the first two values are used.
            # Colour a leaf node if the year is within the target year range;
            # otherwise no colour.
            if ($split_date[0] >= $targets[0] && $split_date[0] <= $targets[1]) {
                $colour_attribute = "[&!color=" . $green . "]"; # green
            } else {
                $colour_attribute = "";
            }
        }
        # append the appropriate colour attribute to the leaf node id
        print $fh "\t", $leaf->id, $colour_attribute, "\n";
    }
    
    print $fh ";\n";
    print $fh "END;\n";
    print $fh "\n";
    close $fh;
    
    # Part 2: write the trees block
    foreach my $node (@nodes) {
        my $boot = $node->bootstrap;
        my $colour_attribute = "";
        my $green = "#00ff00";  # green
        my $red = "#ff0000";    #red
        
        # If the selected feature type is not date, add a colour code to the
        # branch of the node with the matching majority value.
        if ($feature_selector != 1) {
            TREES_MATCH_TARGETS: foreach my $target (@targets) {
                # The majority_percentage tag value has to be copied to
                # another tag called copied_majority_percentage to prevent
                # the original majority_percentage tag value from being
                # modified.
                if ($node->get_tag_values("majority") =~ /$target/i) {
                    $node->set_tag_value("copied_majority_percentage", round($node->get_tag_values("majority_percentage")));
                    $colour_attribute = ",!color=" . $green;    # green
                    last TREES_MATCH_TARGETS;
                } elsif ($target eq "ALL") {
                    # If the default target feature value is used, the
                    # majority feature percentage is still shown, but the
                    # branch is not coloured.
                    $node->set_tag_value("copied_majority_percentage", round($node->get_tag_values("majority_percentage")));
                }
                else {
                    # The majority feature percentage of an unselected branch
                    # is shown as "NA".
                    $node->set_tag_value("copied_majority_percentage", "NA");
                    $colour_attribute = "";
                }
            }
        }
        
        # If the selected feature type is date, add a colour code to the
        # branch of the node with the feature majority value within the
        # target year range. Only the first two values of @targets are used.
        if ($feature_selector == 1) {
            # The majority_percentage tag value has to be copied to another
            # tag called copied_majority_percentage to prevent the original
            # majority_percentage tag value from being modified.
            if ($node->get_tag_values("majority") >= $targets[0] && $node->get_tag_values("majority") <= $targets[1]) {
                $node->set_tag_value("copied_majority_percentage", round($node->get_tag_values("majority_percentage")));
                $colour_attribute = ",!color=" . $green;    # green
            } else {
                $node->set_tag_value("copied_majority_percentage", "NA");
                $colour_attribute = "";
            }
        }
        
        if ($node->is_Leaf) {
            # The majority percentage is stored in the name attribute.
			
			# For a leaf node, add
            # [&!name="majority_percentage",!color=colour_attribute]
			# to the FigTree_parameters tag.
            $node->add_tag_value("FigTree_parameters", "[&!name=\"" . $node->get_tag_values("copied_majority_percentage") . "\"" . $colour_attribute . "]");
            # append the FigTree_parameters tag to its leaf node id
            $node->id($node->id.$node->get_tag_values("FigTree_parameters"));
        } else {
            # For an internal node, colour the node in red if its majority
			# feature value matches the target values and has a bootstrap
			# value greater than or equal to the threshold.
            if ($feature_selector != 1) {
                BOOTSTRAP_MATCH_TARGETS: foreach my $target (@targets) {
                    if ($node->get_tag_values("majority") =~ /$target/i and $boot >= $boot_threshold) {
                        $colour_attribute = ",!color=" . $red;  # red
                        last BOOTSTRAP_MATCH_TARGETS;
                    }
                }
            } else {
                if ($node->get_tag_values("majority") >= $targets[0] and $node->get_tag_values("majority") <= $targets[1] and $boot >= $boot_threshold) {
					# Branches with bootstrap values >= threshold are
					# coloured in red
                    $colour_attribute = ",!color=" . $red;  # red
                }
            }
            # For an internal node, add
            # [&label=bootstrap_value,!name="majority_percentage",
            # !color=colour_attribute] to the FigTree_parameters tag.
            $node->add_tag_value("FigTree_parameters", "&label=" . $boot . ",!name=\"" . $node->get_tag_values("copied_majority_percentage") . "\"" . $colour_attribute);
            
			# Assign the FigTree_parameters tag value to an internal node's
            # bootstrap value slot. This action allows the bootstrap value
            # to be written in square brackets, which are treated as comments
            # in NEXUS format.
            $node->bootstrap($node->get_tag_values("FigTree_parameters"));
        }
    }
    
    # append the output NEXUS tree without the #NEXUS header to the
	# taxa block, and label translation after the taxa block
    my $outtree = Bio::TreeIO->new(-format => "nexus", -file => ">>" . $outfile, -translate => 0, -header => 0);
    $outtree->write_tree($tree_copy);
    
}


# Check the arguments if date is selected to determine majority.
sub _check_date_arguments {
	my ($feature_code, @values) = @_;
	
	if ($feature_code == 1) {
		# Check whether the arguments are integers.
		foreach my $v (@values) {
			if ($v =~ /\D+/) {
				die "Invalid arguments: non-integers (@values) entered to select date.\n";
			}
		}
		# Check whether the start year is greater than the end year.
		if ($values[0] > $values[1]) {
			die "Incorrect order of arguments: start date ($values[0]) is later than end date ($values[1]).\n";
		}
	}
}


# Print a tree summary on the screen.
sub write_tree_summary {
    my ($tree) = @_;
    my @non_binary_nodes = get_non_binary_nodes($tree);
    
    print "Tree summary:\n";
	print "Input file is $ARGV[0]\n";
    print "Total number of nodes is ", scalar $tree->get_nodes, "\n";
    print "Number of leaf nodes is ", scalar $tree->get_leaf_nodes, "\n";
    print "Maximum number of feature fields in node identifiers is $max_fields\n";
    print "Height of tree is ", $tree->height, "\n";
    print "Total branch length is ", $tree->total_branch_length, "\n";
	if (@non_binary_nodes) {
        print "The following nodes are not binary.\n";
        print "node_number", "\t", "id", "\t", "descendant_nodes", "\n";
        foreach my $node (@non_binary_nodes) {
            print $node->get_tag_values("node_number"), "\t", $node->id, "\t";
            my $flag = 0;
            foreach my $descendant ($node->each_Descendent) {
                print "," if ($flag == 1);
                $flag = 1;
                print $descendant->get_tag_values("node_number");
            }
        }
    }
    print "\n";
}

# always return true at the end of the package
1;
