#!/usr/bin/perl -w
#
# tree_extract_OG.pl -- Extracting allelic orthogroups from tree
#
# Author: Yuqian Jiang
# Created: 2025-04-22
# Supported by deepseek R1

use strict;
use warnings;
use Bio::Tree::Compatible;
use Bio::TreeIO;
use Getopt::Long;
use List::Util qw(shuffle);
use Path::Tiny;
use Data::Dumper;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

tree_extract_OG.pl -- Extracting allelic orthogroups from tree

=head1 SYNOPSIS

    tree_extract_OG.pl
        Extracting allelic orthogroups from tree

    Usage:
        perl tree_extract_OG.pl -i <.treefile> -a <taxa.tsv>

    Required:
        -i,--input        STR       input treefile, required
        -a,--add          STR       add taxa info to each node, required (col1:species col2: genes)

    Options:
        -m,--max          FLO       maximum cut-off phylo-tree branch length, default: 0.1
        -t,--taxa         STR       required taxa number cut-off, default: 4
        -r,--reroot       STR       reroot ID if defined, default: unspecified
        -o,--output       STR       output orthogroups in tsv format, default: stdout
        -h,--help                   help information

=cut

GetOptions(
    'input|i=s'     => \(my $input),
    'add|a=s'       => \(my $add_species),
    'output|o=s'    => \(my $output = 'stdout'),
    'max|m=f'       => \(my $max_length = '0.1'),
    'taxa|t=i'      => \(my $required_taxa = '4'),
    'reroot|r=s'    => \(my $reroot_id),
    "h|help"        => sub { Getopt::Long::HelpMessage(0) },
) or Getopt::Long::HelpMessage(1);

if ( !defined $input ) {
    print STDERR "Error: cannot find input tree files.\n";
    die Getopt::Long::HelpMessage(1);
}

if ( !defined $add_species ) {
    print STDERR "Error: cannot find add-info files.\n";
    die Getopt::Long::HelpMessage(1);
}


#----------------------------------------------------------#
# Main program
#----------------------------------------------------------#

our %TAXA;
my $target_node;

# read add species
my @lines = path($add_species) -> lines;
for (@lines) {
    chomp;
    my @arrays = split/\t/, $_;
    $TAXA{$arrays[1]} = $arrays[0];
}

# read input treefile
my $treeio = Bio::TreeIO -> new(
    -file   => $input,
    -format => 'newick'
);
my $tree = $treeio -> next_tree or die "Cannot analyze the tree, wrong newick format.\n";

if (defined $reroot_id) {
    for my $node ($tree -> get_nodes) {
        next unless defined $node -> id;
        if ($node -> id eq $reroot_id) {
            $target_node = $node;
            last;
        }
    }
    die "Cannot find reroot ID '$reroot_id' in all tree labels!\n" if (!defined $target_node);
    eval {
        $tree -> reroot($target_node);
        1;
    } or  die "Reroot failure: $@\n";
}

# pre-calculate each node to root distances
compute_root_distances($tree -> get_root_node);

# collect all gene names
my %all_genes = map { $_ -> id => 1 } $tree -> get_leaf_nodes;

# candidate clusters collection
my @candidate_clusters;
foreach my $node (@{$tree -> Bio::Tree::Compatible::postorder_traversal}) {
    # skip leaf nodes
    next if $node -> is_Leaf;

    # get effective leaf nodes
    my @valid_leaves = get_valid_leaves($node, $max_length);
    next unless @valid_leaves;

    # extract taxa info
    my %taxa_count;
    foreach my $leaf (@valid_leaves) {
        my $taxa = $TAXA{$leaf};
        $taxa_count{$taxa}++;
    }
    # check taxa numbers
    if (scalar keys %taxa_count >= $required_taxa) {
        push @candidate_clusters, {
            genes => \@valid_leaves,
            size  => scalar @valid_leaves,
            taxa  => [keys %taxa_count]
        };
    }
}

# sort according to cluster sizes
@candidate_clusters = sort { $b->{size} <=> $a->{size} } @candidate_clusters;

# filter not intersected clusters
my (%occupied_genes, @final_clusters);
foreach my $cluster (@candidate_clusters) {
    # get unassigned genes
    my @new_genes = grep { !exists $occupied_genes{$_} } @{$cluster->{genes}};
    next unless @new_genes;

    # check taxa numbers
    my %new_taxa;
    foreach my $gene (@new_genes) {
        my $taxon = $TAXA{$gene};
        $new_taxa{$taxon}++;
    }
    next if scalar keys %new_taxa < $required_taxa;

    # record all results
    push @final_clusters, \@new_genes;
    $occupied_genes{$_}++ for @new_genes;
}

# output
my $out_fh;
if ( lc($output) eq "stdout" ) {
    $out_fh = *STDOUT;
}
else {
    open $out_fh, ">", $output;
}

print $out_fh "Gene\tCluster\n";

# mapping clusters and IDs
my $cluster_id = 1;
foreach my $cluster (@final_clusters) {
    print $out_fh "$_\tOG$cluster_id\n" for @$cluster;
    $cluster_id++;
}

# unassigned genes
foreach my $gene (keys %all_genes) {
    next if exists $occupied_genes{$gene};
    print $out_fh "$gene\tnon_OG\n";
}
close $out_fh;


#----------------------------------------------------------#
# Sub program
#----------------------------------------------------------#

sub compute_root_distances {
    my ($node, $parent_dist) = @_;
    # root node parent distance set to 0
    $parent_dist ||= 0;

    # calculate the distance from the node to the root
    my $current_dist = $parent_dist + ($node -> branch_length || 0);
    $node -> {root_dist} = $current_dist;

    # recursive deal with each node
    foreach my $child ($node -> each_Descendent) {
        compute_root_distances($child, $current_dist);
    }
}

sub get_valid_leaves {
    my ($node, $remaining_length) = @_;
    return ($node -> id) if $node -> is_Leaf;

    my @valid;
    foreach my $child ($node -> each_Descendent) {
        my $branch_length = $child -> branch_length || 0;
        next if $branch_length > $remaining_length;

        push @valid, get_valid_leaves($child, $remaining_length - $branch_length);
    }
    return map { ref $_ ? @$_ : $_ } @valid; # flatten the array
}
