#!/usr/bin/perl -w
#
# A script that clustering hhalign all-vs-all results.
#
# Author: Yuqian Jiang
# Created: 2026-08-21
# Updated: Add Getopt::Long, add different cluster methods

use strict;
use warnings;
use Graph::Directed;
use Graph::Undirected;
use Getopt::Long;
use List::Util qw(min max);

=head1 NAME

clan_cluster.pl - clustering domains from hhalign all-vs-all results

=head1 SYNOPSIS

perl clan_cluster.pl -i <hhalign.tsv> -e 1e-10 [options]

Input tsv, one pair per line (header optional):
    Query   Target  E_value Score

Options:
    --input    -i STR   input hhalign tsv
    --evalue   -e STR   E-value cutoff (default: 1e-5)
    --cluster  -c STR   cluster method (default: scc)
                        scc: strongly connected components on directed edges;
                                an edge A->B is kept only if E(A->B) itself passes,
                                so merging needs mutual (direct or transitive) hits;
                                robust to asymmetric scores, no symmetrization
                        cc:  connected components on the undirected graph,
                                pair scores symmetrized first via --strategy
    --strategy -s STR   how to merge E(A->B) and E(B->A), only for -c cc
                        (default: min)
                        min: the better direction passes (permissive, may chain)
                        max: both directions must pass (strict, reciprocal)
                        avg: arithmetic mean
    --help     -h       show help message

Output: Cluster_N <tab> ACC per line to stdout (clusters sorted by size desc),
        summary to stderr.

=cut

GetOptions(
    'help|h'       => sub { Getopt::Long::HelpMessage(0) },
    'input|i=s'    => \( my $input ),
    'evalue|e=s'   => \( my $cutoff_str = '1e-5' ),
    'cluster|c=s'  => \( my $method = 'scc' ),
    'strategy|s=s' => \( my $strategy = 'min' ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $input ) {
    print STDERR "Error: the --input option is required.\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( !-e $input ) {
    die "Error: can't find the input file [$input].\n";
}
if ( $method !~ /^(?:scc|cc)$/ ) {
    die "Error: --cluster must be scc or cc.\n";
}
if ( $strategy !~ /^(?:min|max|avg)$/ ) {
    die "Error: --strategy must be min, max or avg.\n";
}
if ( $method eq 'scc' and $strategy ne 'min' ) {
    print STDERR "Note: --strategy only applies to -c cc, ignored for scc.\n";
}

my $cutoff = 0 + $cutoff_str;

#----------------------------------------------------------#
# Main program
#----------------------------------------------------------#

my %SCORES;        # query -> target -> E-value
my %all_domains;
my $pairs = 0;

open my $IN, '<', $input or die "Cannot open file: $!";
while (<$IN>) {
    chomp;
    next if /^Query/;
    next if /^\s*$/;
    my @arr = split /\t/;
    next if @arr < 3;
    $SCORES{ $arr[0] }{ $arr[1] } = $arr[2];
    $all_domains{ $arr[0] } = 1;
    $all_domains{ $arr[1] } = 1;
    $pairs++;
}
close $IN;

# symmetrized score of an unordered pair
sub sym_score {
    my ( $a, $b ) = @_;
    my $ab = exists $SCORES{$a}{$b} ? $SCORES{$a}{$b} : undef;
    my $ba = exists $SCORES{$b}{$a} ? $SCORES{$b}{$a} : undef;
    if ( defined $ab and defined $ba ) {
        return min( $ab, $ba ) if $strategy eq 'min';
        return max( $ab, $ba ) if $strategy eq 'max';
        return ( $ab + $ba ) / 2;
    }
    elsif ( defined $ab or defined $ba ) {
        # only one direction: max demands reciprocity, others accept the single hit
        return 'inf' if $strategy eq 'max';
        return defined $ab ? $ab : $ba;
    }
    return 'inf';
}

# array of array-refs
my @clusters;

if ( $method eq 'scc' ) {
    my $graph = Graph::Directed->new;
    for my $src ( keys %SCORES ) {
        for my $tgt ( keys %{ $SCORES{$src} } ) {
            next if $src eq $tgt;
            $graph->add_edge( $src, $tgt ) if $SCORES{$src}{$tgt} + 0 < $cutoff;
        }
    }
    push @clusters, $_ for $graph->strongly_connected_components;

    # isolated nodes as singleton clusters
    my %in_graph = map { $_ => 1 } $graph->vertices;
    for my $node ( keys %all_domains ) {
        push @clusters, [$node] if !$in_graph{$node};
    }
}
else {
    my $graph = Graph::Undirected->new;
    $graph->add_vertex($_) for keys %all_domains;
    for my $src ( keys %SCORES ) {
        for my $tgt ( keys %{ $SCORES{$src} } ) {
            next if $src eq $tgt;
            next if $src gt $tgt;    # visit each unordered pair once
            $graph->add_edge( $src, $tgt ) if sym_score( $src, $tgt ) < $cutoff;
        }
    }
    push @clusters, $_ for $graph->connected_components;
}

# sort clusters by size desc and renumber
@clusters = sort { scalar(@$b) <=> scalar(@$a) } @clusters;

for my $i ( 0 .. $#clusters ) {
    my $num = $i + 1;
    print "Cluster_$num\t$_\n" for @{ $clusters[$i] };
}

my $singletons = grep { scalar(@$_) == 1 } @clusters;
printf STDERR
    $method, $method eq 'scc' ? '-' : $strategy, $cutoff_str,
    scalar(keys %all_domains), $pairs, scalar @clusters,
    scalar @{ $clusters[0] }, $singletons;

__END__
