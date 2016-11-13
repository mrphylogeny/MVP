# MVP.pm
Perl package that determines the majority feature values of the nodes of a phylogenetic tree by simple majority rule and parsimony.
The pacakage requires a Newick tree as input. The node identifiers of the Newick tree should contain values in @-delimited fields. The values at a particular field represent a feature of that tree. The node identifiers have to be specifically formatted as accession_number@date@feature1@feature2@feature3@..., in which the accession number and date are mandatory, and follows with an arbitrary number of optional features.
