#!/usr/bin/perl -w
#
# graph_cluster.pl -- Graph clustering on edge lists
#
# Author: Yuqian Jiang
# Created: 2026-08-21
# Updated: 2026-08-21 Added SCC-based clustering
# Updated: 2026-08-21 Generalized: edge-list input, multiple cluster methods (cc/wcc/scc/lpa/louvain/mcl)

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max sum);

Getopt::Long::Configure('no_ignore_case');

my ( $input, $out, $dir, $method, $strategy, $cutoff_str, $greater,
    $inflation, $resolution, $seed, $header, $verbose, $help );

GetOptions(
    'input|i=s'        => \$input,
    'out|o=s'          => \$out,
    'dir|d=s'          => \$dir,
    'cluster|c=s'      => \$method,
    'strategy|s=s'     => \$strategy,
    'evalue|e=s'       => \$cutoff_str,
    'greater|g'        => \$greater,
    'inflation|r=f'    => \$inflation,
    'resolution|R=f'   => \$resolution,
    'seed=i'           => \$seed,
    'header|H'         => \$header,
    'verbose|v'        => \$verbose,
    'help|h'           => \$help,
) or die "Use -h for help.\n";

my $usage = <<'__GUIDE__';
Graph clustering on an edge list (pure Perl, no external dependency).

Input tsv, one edge per line ('#' comments and extra columns ignored):
    col1: node1   col2: node2   col3: value (optional)
    col3 is an E-value-like weight (smaller = better) by default;
    a two-column line is treated as unweighted (value = 1).

Usage:
    perl graph_cluster.pl -i <edges.tsv> -c <method> [options]
    perl graph_cluster.pl -i <edges.tsv> -e 1e-10 -c scc
    perl graph_cluster.pl -i <mmseqs_cluster.tsv> -c cc

Required:
    -i,--input     <file>  edge list, '-' for STDIN

Cluster methods (-c, default: cc):
    cc      undirected connected components (single linkage),
            pair scores symmetrized first via -s
    wcc     weakly connected components, direction ignored entirely
            (an edge exists if ANY direction passes the cutoff)
    scc     strongly connected components on directed edges;
            merging needs mutual (direct or transitive) evidence,
            robust to asymmetric scores, -s not used
    lpa     label propagation (fast, stochastic, --seed controls it)
    louvain modularity optimization, communities of dense links,
            -R resolution (default 1.0, larger gives more clusters)
    mcl     Markov clustering, -r inflation (default 2.0,
            larger gives more, tighter clusters)

Edge filtering and symmetrizing:
    -e,--evalue   <num>   cutoff, default: keep all edges;
                          pass if value < cutoff (or >= with -g)
    -g,--greater          value is similarity-like (greater = better)
    -s,--strategy <str>   how to merge the two directions of a pair
                          (undirected methods only, default: min)
                          min | max | mean
                          for E-values, min = any direction passes,
                          max = both directions must pass (reciprocal);
                          for similarity (-g) the two are swapped

Output (results never mixed with running messages):
    -o,--out      <file>  Cluster_N <tab> node, one per line
                          (default: STDOUT; messages go to STDERR)
    -d,--dir      <dir>   additionally write cluster_N.lst per cluster
                          and single.all.lst for singletons
                          (compatible with cluster_undirect.pl)

Others:
    --seed <int>  random seed for lpa/louvain (default: 42)
    -H,--header   skip the first line of input
    -v,--verbose  print parsing details to STDERR
    -h,--help     show this help

__GUIDE__

die $usage if $help;

#----------------------------------------------------------#
# Init and checks
#----------------------------------------------------------#

$method    = 'cc'        unless defined $method;
$strategy  = 'min'       unless defined $strategy;
$inflation = 2.0         unless defined $inflation;
$resolution= 1.0         unless defined $resolution;
$seed      = 42          unless defined $seed;

if ( !defined $input ) {
    print STDERR "Error: the --input option is required.\n";
    die $usage;
}
my $fh_in;
if ( $input eq '-' ) { $fh_in = \*STDIN; }
elsif ( -e $input )  { open $fh_in, '<', $input or die "Error: can't open [$input]\n"; }
else { die "Error: can't find the input file [$input]\n"; }

if ( $method !~ /^(?:cc|wcc|scc|lpa|louvain|mcl)$/ ) {
    die "Error: --cluster must be cc, wcc, scc, lpa, louvain or mcl.\n";
}
if ( $strategy !~ /^(?:min|max|mean)$/ ) {
    die "Error: --strategy must be min, max or mean.\n";
}
if ( $method eq 'scc' and $strategy ne 'min' ) {
    print STDERR "Note: --strategy is not used by -c scc (directed edges).\n";
}

my $cutoff = defined $cutoff_str ? 0 + $cutoff_str : undef;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#

sub pass_cutoff {
    my ($v) = @_;
    return 1 if !defined $cutoff;
    return $v >= $cutoff ? 1 : 0 if $greater;
    return $v < $cutoff ? 1 : 0;
}

sub sym_of {
    my ($vals) = @_;    # arrayref of all observed values of an unordered pair
    return min(@$vals) if $strategy eq 'min';
    return max(@$vals) if $strategy eq 'max';
    return sum(@$vals) / scalar(@$vals);
}

# deterministic seeded shuffle
sub shuffle_seeded {
    my ( $nodes, $s ) = @_;
    my @arr = @$nodes;
    srand($s);
    for my $i ( reverse 1 .. $#arr ) {
        my $j = int( rand( $i + 1 ) );
        @arr[ $i, $j ] = @arr[ $j, $i ];
    }
    return @arr;
}

# union-find connected components over [a, b] edge list
sub union_find {
    my ($edges) = @_;
    my %parent;
    my $find = sub {
        my $x = shift;
        my @path;
        while ( $parent{$x} ne $x ) { push @path, $x; $x = $parent{$x}; }
        $parent{$_} = $x for @path;
        return $x;
    };
    for my $e (@$edges) {
        my ( $a, $b ) = @$e;
        $parent{$a} = $a unless exists $parent{$a};
        $parent{$b} = $b unless exists $parent{$b};
        my ( $ra, $rb ) = ( $find->($a), $find->($b) );
        $parent{$rb} = $ra if $ra ne $rb;
    }
    my %groups;
    for my $n ( keys %parent ) { push @{ $groups{ $find->($n) } }, $n; }
    return values %groups;
}

# iterative DFS collecting postorder (finish order)
sub postorder {
    my ( $adj, $start, $seen, $order ) = @_;
    return if $seen->{$start};
    my @stack = [ $start, [ sort keys %{ $adj->{$start} || {} } ] ];
    $seen->{$start} = 1;
    while (@stack) {
        my $frame = $stack[-1];
        if ( @{ $frame->[1] } ) {
            my $nxt = shift @{ $frame->[1] };
            next if $seen->{$nxt};
            $seen->{$nxt} = 1;
            push @stack, [ $nxt, [ sort keys %{ $adj->{$nxt} || {} } ] ];
        }
        else {
            push @$order, $frame->[0];
            pop @stack;
        }
    }
}

# Kosaraju strongly connected components
sub kosaraju {
    my ($adj) = @_;
    my ( %radj, %seen, @order );
    for my $a ( keys %$adj ) {
        $radj{$_}{$a} = 1 for keys %{ $adj->{$a} };
    }
    postorder( $adj, $_, \%seen, \@order ) for sort keys %$adj;
    my @comp;
    %seen = ();
    for my $n ( reverse @order ) {
        next if $seen{$n};
        my @collect;
        postorder( \%radj, $n, \%seen, \@collect );
        push @comp, [@collect];
    }
    return @comp;
}

# label propagation on weighted undirected adjacency {a}{b}=w
sub lpa {
    my ($adj) = @_;
    my %label = map { $_ => $_ } keys %$adj;
    my $changed = 1;
    my $iter = 0;
    while ($changed) {
        $changed = 0;
        $iter++;
        for my $n ( shuffle_seeded( [ sort keys %$adj ], $seed + $iter ) ) {
            my %vote;
            $vote{ $label{$_} } += $adj->{$n}{$_} for keys %{ $adj->{$n} };
            next if !%vote;
            my ($best) = sort { $vote{$b} <=> $vote{$a} || $a cmp $b } keys %vote;
            if ( $best ne $label{$n} ) { $label{$n} = $best; $changed = 1; }
        }
        last if $iter >= 100;
    }
    my %groups;
    push @{ $groups{ $label{$_} } }, $_ for keys %$adj;
    return values %groups;
}

# Louvain modularity optimization on weighted undirected adjacency {a}{b}=w
sub louvain {
    my ($adj) = @_;

    # level graph, symmetric, self-loops stored once in g{a}{a};
    # degree k{a} = sum_b g{a}{b} + g{a}{a}, m = sum(k)/2
    my %g;
    for my $a ( keys %$adj ) {
        $g{$a}{$_} = $adj->{$a}{$_} for keys %{ $adj->{$a} };
    }
    my %where = map { $_ => $_ } keys %$adj;    # orig node -> current supernode

    my $level = 0;
    while ( $level++ < 100 ) {

        # degrees
        my ( %k, $sum_k );
        for my $a ( keys %g ) {
            my $s = 0;
            $s += $_ for values %{ $g{$a} };
            $s += $g{$a}{$a} if exists $g{$a}{$a};
            $k{$a} = $s;
            $sum_k += $s;
        }
        my $m = $sum_k / 2;
        last if $m <= 0;

        # one community per node
        my %comm      = map { $_ => $_ } keys %g;
        my %sigma_tot;
        $sigma_tot{ $comm{$_} } += $k{$_} for keys %g;

        # local moving
        my ( $moved, $iter ) = ( 1, 0 );
        while ($moved) {
            $moved = 0;
            $iter++;
            for my $a ( shuffle_seeded( [ sort keys %g ], $seed + $iter + $level ) ) {
                my %wc;    # weights to neighbor communities (self excluded)
                for my $b ( keys %{ $g{$a} } ) {
                    next if $b eq $a;
                    $wc{ $comm{$b} } += $g{$a}{$b};
                }
                my $ca = $comm{$a};
                $sigma_tot{$ca} = ( $sigma_tot{$ca} // 0 ) - $k{$a};    # remove a

                my ( $best, $best_gain ) = ( $ca, -'inf' );
                my %cands = %wc;
                $cands{$ca} = $wc{$ca} // 0;
                for my $c ( keys %cands ) {
                    my $gain = $cands{$c} - $resolution * ( $sigma_tot{$c} // 0 ) * $k{$a} / $m;
                    if ( $gain > $best_gain ) { $best_gain = $gain; $best = $c; }
                }
                if ( $best_gain < 0 ) {    # isolate into a fresh community
                    $best = "__solo_$a";
                    $sigma_tot{$best} = 0;
                }
                $comm{$a} = $best;
                $sigma_tot{$best} += $k{$a};
                $moved++ if $best ne $ca;
            }
            last if $iter >= 50;
        }

        # aggregate
        my %members;
        push @{ $members{ $comm{$_} } }, $_ for keys %g;
        last if scalar(keys %members) == scalar( keys %g );   # converged

        $where{$_} = $comm{ $where{$_} } for keys %where;

        my %ng;
        for my $a ( sort keys %g ) {
            for my $b ( sort keys %{ $g{$a} } ) {
                next if $a gt $b;    # each unordered edge once
                my ( $ca, $cb ) = ( $comm{$a}, $comm{$b} );
                my $w = $g{$a}{$b};
                if ( $ca eq $cb ) { $ng{$ca}{$ca} += $w }
                else { $ng{$ca}{$cb} += $w; $ng{$cb}{$ca} += $w; }
            }
        }
        %g = %ng;
    }

    my %groups;
    push @{ $groups{ $where{$_} } }, $_ for keys %where;
    return values %groups;
}

# Markov clustering on weighted undirected adjacency, inflation $r
sub mcl {
    my ( $adj, $r ) = @_;
    my @nodes = sort keys %$adj;

    # normalized weights in (0,1], self-loop = mean weight
    my $maxw = 0;
    for my $a (@nodes) {
        for my $b ( keys %{ $adj->{$a} } ) { $maxw = $adj->{$a}{$b} if $adj->{$a}{$b} > $maxw; }
    }
    my $n_edges = 0;
    $n_edges += scalar keys %{ $adj->{$_} } for @nodes;
    my $sum_w = 0;
    $sum_w += $_ for map { values %{ $adj->{$_} } } @nodes;
    my $selfw = $n_edges ? $sum_w / $n_edges / $maxw : 1;

    # row-stochastic matrix M{a}{b}
    my %M;
    for my $a (@nodes) {
        $M{$a}{$a} = $selfw;
        $M{$a}{$_} = $adj->{$a}{$_} / $maxw for keys %{ $adj->{$a} };
        my $rowsum = sum( values %{ $M{$a} } );
        $M{$a}{$_} /= $rowsum for keys %{ $M{$a} };
    }

    my $converged = 0;
    for my $round ( 1 .. 100 ) {

        # expand: E = M x M
        my %E;
        for my $a (@nodes) {
            for my $b ( keys %{ $M{$a} } ) {
                my $pab = $M{$a}{$b};
                next if $pab < 1e-9;
                for my $c ( keys %{ $M{$b} } ) {
                    $E{$a}{$c} += $pab * $M{$b}{$c};
                }
            }
        }

        # inflate, prune, renormalize, check convergence
        my $maxdiff = 0;
        my %oldm;
        for my $a (@nodes) { $oldm{$a} = { %{ $M{$a} || {} } } }
        %M = ();
        for my $a (@nodes) {
            my $rowsum = 0;
            for my $c ( keys %{ $E{$a} } ) {
                my $v = $E{$a}{$c}**$r;
                next if $v < 1e-8;
                $M{$a}{$c} = $v;
                $rowsum += $v;
            }
            next unless $rowsum > 0;
            $M{$a}{$_} /= $rowsum for keys %{ $M{$a} };
        }
        for my $a (@nodes) {
            my %u = map { $_ => 1 } ( keys %{ $M{$a} }, keys %{ $oldm{$a} } );
            for my $c ( keys %u ) {
                my $d = abs( ( $M{$a}{$c} // 0 ) - ( $oldm{$a}{$c} // 0 ) );
                $maxdiff = $d if $d > $maxdiff;
            }
        }
        if ( $maxdiff < 1e-5 ) { $converged = $round; last; }
    }
    print STDERR "MCL " . ( $converged ? "converged in $converged rounds" : "hit round cap (100)" ) . "\n"
        if $verbose;

    # clusters: group nodes by their argmax attractor
    my %attract;
    for my $n (@nodes) {
        my ($best)
            = sort { ( $M{$n}{$b} <=> $M{$n}{$a} ) || $a cmp $b } keys %{ $M{$n} };
        $attract{$n} = $best;
    }
    my %groups;
    push @{ $groups{ $attract{$_} } }, $_ for @nodes;
    return values %groups;
}

#----------------------------------------------------------#
# Main program: parse edges
#----------------------------------------------------------#

my ( %DIR, %NODE );
my ( $line_n, $warned ) = ( 0, 0 );

LINE: while (<$fh_in>) {
    chomp;
    next if /^\s*$/;
    next if /^#/;
    $line_n++;
    if ( $header && $line_n == 1 ) { $header = 0; next; }
    my @f = split /\t/;
    next LINE if @f < 2;
    my ( $a, $b ) = ( $f[0], $f[1] );
    my $v = 1;
    if ( defined $f[2] && $f[2] ne '' ) {
        if ( $f[2] =~ /^-?\d+(\.\d+)?([eE][+-]?\d+)?$/ ) { $v = $f[2] + 0 }
        elsif ( !$warned ) {
            print STDERR "Warning: non-numeric col3 [$f[2]] at line $., treated as unweighted.\n";
            $warned = 1;
        }
    }
    $NODE{$a} = 1;
    $NODE{$b} = 1;
    next if $a eq $b;    # self-loops only register the node

    # directed store, keep the best duplicate per direction
    if ( !exists $DIR{$a}{$b} ) { $DIR{$a}{$b} = $v }
    else {
        $DIR{$a}{$b} = min( $DIR{$a}{$b}, $v ) unless $greater;
        $DIR{$a}{$b} = max( $DIR{$a}{$b}, $v ) if $greater;
    }
}

#----------------------------------------------------------#
# Main program: filter and build graphs
#----------------------------------------------------------#

# all observed values per unordered pair (pre-filter, both directions)
my ( %UND, %NDIR );    # values, and which directions were observed
for my $src ( keys %DIR ) {
    for my $tgt ( keys %{ $DIR{$src} } ) {
        my $key = $src le $tgt ? "$src\t$tgt" : "$tgt\t$src";
        push @{ $UND{$key} }, $DIR{$src}{$tgt};
        $NDIR{$key}{$src} = 1;
    }
}

# reciprocal mode: the strategy demands BOTH directions to pass
my $reciprocal
    = ( !$greater && $strategy eq 'max' ) || ( $greater && $strategy eq 'min' )
        ? 1 : 0;

my ( @dedges, @uedges, %USYM );    # directed / undirected edges, sym values
if ( $method eq 'scc' ) {
    for my $a ( keys %DIR ) {
        for my $b ( keys %{ $DIR{$a} } ) {
            push @dedges, [ $a, $b ] if pass_cutoff( $DIR{$a}{$b} );
        }
    }
}
else {
    for my $key ( keys %UND ) {
        my ( $a, $b ) = split /\t/, $key;
        my $vals = $UND{$key};
        my $ok;
        if ( $method eq 'wcc' ) { $ok = ( grep { pass_cutoff($_) } @$vals ) ? 1 : 0 }
        else {
            $ok = pass_cutoff( sym_of($vals) );
            $ok = 0 if $ok and $reciprocal and keys %{ $NDIR{$key} } < 2;
        }
        next unless $ok;
        push @uedges, [ $a, $b ];
        $USYM{$key} = sym_of($vals);
    }
}

print STDERR "Parsed $line_n lines; ",
    ( $method eq 'scc' ? scalar(@dedges) . " directed edges" : scalar(@uedges) . " undirected pairs" ),
    " pass the cutoff.\n"
    if $verbose;

#----------------------------------------------------------#
# Main program: cluster
#----------------------------------------------------------#

my @clusters;
my %linked;    # nodes covered by at least one passing edge

if ( $method eq 'scc' ) {
    my %adj_h;
    push @{ $adj_h{ $_->[0] } }, $_->[1] for @dedges;
    my %adj = map { $_ => { map { $_ => 1 } @{ $adj_h{$_} } } } keys %adj_h;
    @clusters = kosaraju( \%adj );
    $linked{$_} = 1 for map { @$_ } @dedges;
}
else {
    # weighted undirected adjacency
    my %adj;
    for my $e (@uedges) {
        my ( $a, $b ) = @$e;
        my $key = $a le $b ? "$a\t$b" : "$b\t$a";
        $linked{$a} = 1;
        $linked{$b} = 1;
        if ( $method eq 'cc' || $method eq 'wcc' ) {
            $adj{$a}{$b} = 1;    # existence only
            $adj{$b}{$a} = 1;
        }
        else {
            my $w = $greater ? $USYM{$key} : -log( $USYM{$key} > 0 ? $USYM{$key} : 1e-300 ) / log(10);
            $adj{$a}{$b} = $w;
            $adj{$b}{$a} = $w;
        }
    }
    if    ( $method eq 'cc' || $method eq 'wcc' ) { @clusters = union_find( \@uedges ) }
    elsif ( $method eq 'lpa' )                    { @clusters = lpa( \%adj ) }
    elsif ( $method eq 'louvain' )                { @clusters = louvain( \%adj ) }
    elsif ( $method eq 'mcl' )                    { @clusters = mcl( \%adj, $inflation ) }
}

# singletons for nodes without any passing edge
for my $n ( sort keys %NODE ) {
    push @clusters, [$n] if !$linked{$n};
}

@clusters = sort { scalar(@$b) <=> scalar(@$a) || $a->[0] cmp $b->[0] } @clusters;

#----------------------------------------------------------#
# Main program: output
#----------------------------------------------------------#

my $to_file = defined $out && $out ne 'stdout' ? 1 : 0;
if ($to_file) {
    open my $fh, '>', $out or die "Error: can't write [$out]\n";
    for my $i ( 0 .. $#clusters ) {
        print {$fh} "Cluster_", $i + 1, "\t$_\n" for @{ $clusters[$i] };
    }
    close $fh;
}
else {
    for my $i ( 0 .. $#clusters ) {
        print "Cluster_", $i + 1, "\t$_\n" for @{ $clusters[$i] };
    }
}

# per-cluster lst files, cluster_undirect.pl style
if ( defined $dir ) {
    require File::Path;
    File::Path::make_path($dir) if !-d $dir;
    my ( $i, @single );
    for my $c (@clusters) {
        if ( @$c == 1 ) { push @single, @$c }
        else {
            $i++;
            open my $LST, '>', "$dir/cluster_$i.lst" or die "Error: can't write to [$dir]\n";
            print {$LST} "$_\n" for @$c;
            close $LST;
        }
    }
    open my $SGL, '>', "$dir/single.all.lst" or die "Error: can't write to [$dir]\n";
    print {$SGL} "$_\n" for @single;
    close $SGL;
}

my $singletons = grep { scalar(@$_) == 1 } @clusters;
printf STDERR
    "Method: %s%s, cutoff: %s, %d nodes, %d clusters (largest: %d, singletons: %d).\n",
    $method,
    $method =~ /^(?:cc|lpa|louvain|mcl)$/ ? " (strategy: $strategy)" : '',
    defined $cutoff_str ? $cutoff_str : 'none',
    scalar( keys %NODE ), scalar @clusters,
    scalar @{ $clusters[0] }, $singletons;

__END__
