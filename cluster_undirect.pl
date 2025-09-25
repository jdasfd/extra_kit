#!/usr/bin/perl -w
#
# cluster_undirect.pl -- Using undirect graph to get clusters
#
# Author: Yuqian Jiang
# Created: 2025-04-24

use strict;
use warnings;
use Path::Tiny;
use Graph::Undirected;
use Getopt::Long;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

cluster_undirect.pl -- Using undirect graph to get clusters

=head1 SYNOPSIS

    cluster_undirect.pl
        Using undirect graph to get clusters

    Usage:
        perl cluster_undirect.pl <cluster.tsv>

    Required:
        -i,--input        STR       input cluster.tsv from mmseqs2 (col1: gene1 col2: gene2)
        -d,--dir          STR       the path dir for all results saved

=cut

GetOptions(
    'input|i=s'     => \(my $input),
    'd|dir=s'       => \(my $dir),
    "h|help"        => sub { Getopt::Long::HelpMessage(0) }
) or Getopt::Long::HelpMessage(1);

if ( !defined $input ) {
    print STDERR "Error: cannot find input files.\n";
    die Getopt::Long::HelpMessage(1);
}
if ( !defined $dir ) {
    print STDERR "Error: undefined dir.\n";
    die Getopt::Long::HelpMessage(1);
}

our $g = Graph::Undirected -> new;
our $i = 1;
my @single;

#----------------------------------------------------------#
# Main program
#----------------------------------------------------------#

# read_in files
my @lines = path($input) -> lines;

# split clusters
for (@lines) {
    chomp;
    my @gene = split/\t/, $_;
    $gene[0] eq $gene[1] ? $g -> add_vertex($gene[0]) : $g -> add_edge($gene[0], $gene[1]);
}

# connect all nodes
my @ccs = $g -> connected_components;
# sort connected components (ccs)
@ccs = map { $_->[0] }
    sort { $b->[1] <=> $a->[1] }
        map { [ $_, scalar( @{$_} ) ] } @ccs;

for my $cc (@ccs) {
    my $count = scalar @{$cc};
    if ($count == 1) {
        push @single, @{$cc};
    }
    else {
        path(qq{$dir/cluster_$i.lst}) -> spew(map {qq{$_\n}} @{$cc});
        $i++;
    }
}

path(qq{$dir/single.all.lst}) -> spew(map {qq{$_\n}} @single);
