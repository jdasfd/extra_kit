#!/usr/bin/perl -w
#
# A script that extracting stockholm alignments (Pfam-A.seed)
# or hmm profiles (Pfam-A.hmm) by a list of accessions.
#
# Author: Yuqian Jiang
# Created: 2025-07-20
# Updated: 2026-08-20 Add GetOptions, multi-process extraction and index-based random access.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use File::Path qw(make_path);

=head1 NAME

sto_extract.pl - extracting stockholm alignments or hmm profiles from Pfam files

=head1 SYNOPSIS

perl sto_extract.pl -l <acc.lst> -s <Pfam-A.seed|Pfam-A.hmm> [options]

Options:
    --list   -l STR   a list file of accessions, one per line
                        (e.g. PF00931.29 or PF00931, both supported)
    --seed   -s STR   the Pfam-A.seed (sto blocks) or Pfam-A.hmm file
    --dir    -d STR   output dir, each acc will be written as
                        <dir>/<ACC>.sto or <ACC>.hmm; without it, cat to STDOUT
    --thread -t INT   number of processes for per-file output (default: 1),
                        only works with --dir
    --index  -x STR   index file for random access (default: <seed>.idx),
                        will be used automatically if it exists
    --mkidx          build the index from the seed file and exit
    --help   -h      show help message

Examples:
    # build an index once, speeds up all later runs (low memory)
    perl sto_extract.pl -s Pfam-A.seed --mkidx

    # extract and cat all seed alignments to one file
    perl sto_extract.pl -l CL0023.hmm.lst -s Pfam-A.seed > CL0023.sto

    # extract each seed alignment into a dir with 10 processes
    perl sto_extract.pl -l CL0023.hmm.lst -s Pfam-A.seed -d sto/ -t 10

=cut

GetOptions(
    'help|h'     => sub { Getopt::Long::HelpMessage(0) },
    'list|l=s'   => \( my $list ),
    'seed|s=s'   => \( my $seed ),
    'dir|d=s'    => \( my $dir ),
    'thread|t=i' => \( my $thread = 1 ),
    'index|x=s'  => \( my $index ),
    'mkidx'      => \( my $mkidx ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $seed or !-e $seed ) {
    die "Error: can't find the seed file [$seed].\n";
}
if ( !defined $list and !$mkidx ) {
    print STDERR "Error: one of --list or --mkidx is required.\n";
    die Getopt::Long::HelpMessage(1);
}
if ( $thread > 1 and !defined $dir ) {
    print STDERR "Warning: --thread only works with --dir, force to 1.\n";
    $thread = 1;
}

my $index_file = defined $index ? $index : "$seed.idx";
my $ext = $seed =~ m/\.hmm/ ? 'hmm' : 'sto';

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#

# key line: '#=GF AC PF00931.29' in seed, or 'ACC PF00931.29' in hmm
sub acc_of {
    my ($line) = @_;
    return $1 if $line =~ /^(?:#=GF\s+AC|^ACC)\s+(\S+)/;
    return;
}

sub build_index {
    my ( $seed, $index_file ) = @_;
    open my $IN,  '<', $seed;
    open my $OUT, '>', $index_file;
    my ( $acc, $start ) = ( undef, 0 );
    while (<$IN>) {
        if (/^\/\//) {
            my $end = tell($IN);
            print {$OUT} "$acc\t$start\t", $end - $start, "\n" if defined $acc;
            $start = $end;
            $acc   = undef;
        }
        elsif ( defined( my $hit = acc_of($_) ) ) {
            $acc = $hit;
        }
    }
    close $OUT;
    close $IN;
    print STDERR "Index built: $index_file\n";
    return;
}

sub strip_ver {
    my ($acc) = @_;
    ( my $unver = $acc ) =~ s/\.\d+$//;
    return $unver;
}

#----------------------------------------------------------#
# Main program
#----------------------------------------------------------#

# build the index and exit
if ($mkidx) {
    build_index( $seed, $index_file );
    exit;
}

# read the list file
my @acc_lst;
open my $LST_IN, '<', $list;
while (<$LST_IN>) {
    chomp;
    next if /^\s*$/;
    push @acc_lst, $_;
}
close $LST_IN;

# load data: index (random access) or one scan keeping wanted blocks only
our %SPACE;      # versioned ACC -> content (scan mode) or [offset, length] (index mode)
our %ALIAS;      # unversioned ACC -> versioned ACC
my $use_index = -e $index_file;

if ($use_index) {
    print STDERR "Using index: $index_file\n";
    open my $IDX_IN, '<', $index_file;
    while (<$IDX_IN>) {
        chomp;
        my ( $acc, $off, $len ) = split /\t/;
        $SPACE{$acc} = [ $off, $len ];
        $ALIAS{ strip_ver($acc) } = $acc;
    }
    close $IDX_IN;
}
else {
    my %want;
    for my $acc (@acc_lst) {
        $want{$acc} = 1;
        $want{ strip_ver($acc) } = 1;
    }
    my ( $acc, @buf );
    open my $SEED_IN, '<', $seed;
    while (<$SEED_IN>) {
        if (/^\/\//) {
            if ( defined $acc and ( $want{$acc} or $want{ strip_ver($acc) } ) ) {
                my $content = join( '', @buf ) . "//\n";
                $SPACE{$acc} = $content;
                $ALIAS{ strip_ver($acc) } = $acc;
            }
            ( $acc, @buf ) = ();
        }
        else {
            push @buf, $_;
            my $hit = acc_of($_);
            $acc = $hit if defined $hit;
        }
    }
    close $SEED_IN;
}

# resolve a list item to the versioned ACC
sub resolve {
    my ($item) = @_;
    return $item if exists $SPACE{$item};
    return $ALIAS{ strip_ver($item) } if exists $ALIAS{ strip_ver($item) };
    return;
}

# fetch the block content of a resolved ACC
sub fetch {
    my ( $real, $fh ) = @_;
    if ($use_index) {
        my ( $off, $len ) = @{ $SPACE{$real} };
        seek $fh, $off, 0;
        my $buf;
        read $fh, $buf, $len;
        return $buf;
    }
    else {
        return $SPACE{$real};
    }
}

make_path($dir) if defined $dir and !-d $dir;

# resolve all items first, report missing ones
my @jobs;    # [list_item, resolved_acc]
my @missing;
for my $item (@acc_lst) {
    my $real = resolve($item);
    if ( !defined $real ) {
        print STDERR "Warning: can't find [$item], skipped.\n";
        push @missing, $item;
    }
    else {
        push @jobs, [ $item, $real ];
    }
}

my $ok_count = 0;

if ( defined $dir and $thread > 1 ) {

    # multi-process mode, round-robin over jobs
    my @pids;
    for my $worker ( 0 .. $thread - 1 ) {
        my $pid = fork();
        die "Error: fork failed.\n" if !defined $pid;
        if ( $pid == 0 ) {
            open my $FH, '<', $seed or die "Error: can't open [$seed].\n";
            for my $i ( grep { $_ % $thread == $worker } 0 .. $#jobs ) {
                my $real = $jobs[$i][1];
                open my $OUT, '>', "$dir/$real.$ext";
                print {$OUT} fetch( $real, $FH );
                close $OUT;
            }
            close $FH;
            exit 0;
        }
        push @pids, $pid;
    }
    waitpid( $_, 0 ) for @pids;

    # recount by file existence
    for my $job (@jobs) {
        $ok_count++ if -e "$dir/$job->[1].$ext";
    }
}
else {

    # sequential mode
    my $FH;
    open $FH, '<', $seed or die "Error: can't open [$seed].\n" if $use_index;
    for my $job (@jobs) {
        my ( $item, $real ) = @$job;
        my $content = fetch( $real, $FH );
        if ( defined $dir ) {
            open my $OUT, '>', "$dir/$real.$ext";
            print {$OUT} $content;
            close $OUT;
        }
        else {
            print $content;
        }
        $ok_count++;
    }
    close $FH if $use_index;
}

print STDERR "Total $ok_count of ", scalar(@acc_lst), " items extracted";
print STDERR ", ", scalar(@missing), " missing" if @missing;
print STDERR ".\n";

__END__
