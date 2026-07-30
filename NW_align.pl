#!/usr/bin/perl -w
#
# NW_align.pl -- Needleman-Wunsch global alignment algorithm for DNA / Protein sequences
#
# Author: Yuqian Jiang
# Created: 2026-07-31
# Supported: GLM 5.2
# Referenced:
#   - Needleman_Wunsch.pl (https://github.com/GouXiangJian/two_seq_alignment/blob/master/Needleman_Wunsch.pl)
# Change logs
# 2026-07-31: merge strengths of both references; add DNA/Protein dual mode and BLOSUM62 substitution matrix support

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

my ($fasta, $file1, $file2, $seq1_in, $seq2_in, $type, $matrix,
    $output, $match, $unmatch, $gap, $detailed, $verbose, $help);

GetOptions(
    'fasta|f=s'    => \$fasta,
    'file1|f1=s'   => \$file1,
    'file2|f2=s'   => \$file2,
    'seq1|s1=s'    => \$seq1_in,
    'seq2|s2=s'    => \$seq2_in,
    'type|t=s'     => \$type,
    'matrix|mtx=s' => \$matrix,
    'match|m=i'    => \$match,
    'unmatch|u=i'  => \$unmatch,
    'gap|g=i'      => \$gap,
    'output|o=s'   => \$output,
    'detailed|d'   => \$detailed,
    'verbose|v'    => \$verbose,
    'help|h+'      => \$help,
) or die "Use -h for help.\n";

my $usage = <<'__GUIDE__';
Global alignment between two sequences using Needleman-Wunsch algorithm.

Usage:
    perl NW_align.pl -f <fasta>                       [options]
    perl NW_align.pl -s1 <str> -s2 <str>              [options]
    perl NW_align.pl -f1 <f1> -f2 <f2>                [options]

Required Input (one mode only):
    -f  <file> : FASTA file containing exactly 2 sequences
    -s1 <str>  : raw sequence 1   (use together with -s2)
    -s2 <str>  : raw sequence 2   (use together with -s1)
    -f1 <file> : plain text file with sequence 1
    -f2 <file> : plain text file with sequence 2

Sequence type:
    -t|--type   <str> : dna | protein | auto  (default: auto, DNA if only ACGTUN chars, else protein)

Scoring:
    -mtx|--matrix <str> : simple | blosum62  (default: simple)
                        protein + blosum62 uses substitution matrix
    -m|--match   <int>  : match    score (simple mode, default: 1 DNA / 5 protein)
    -u|--unmatch <int>  : mismatch score (simple mode, default: -1 DNA / -4 protein)
    -g|--gap     <int>  : gap      penalty (default: -2 DNA / -8 protein)

Output:
    -o|--output  <file> : output file (default: STDOUT)
    -d|--detailed       : detailed report (BLAST-like format)
    -v|--verbose        : also print score matrix and path (with -d)
    -h|--help           : show this help

__GUIDE__

die $usage if $help;

# Defaults for type / matrix
$type   = 'auto'   unless defined $type;
$matrix = 'simple' unless defined $matrix;
$type   = lc $type;
$matrix = lc $matrix;

# Resolve sequences based on the chosen input mode
my ($seq1, $id1, $seq2, $id2);
if (defined $fasta) {
    die "Cannot mix -f with -s1/-s2/-f1/-f2.\n"
        if defined $seq1_in || defined $seq2_in || defined $file1 || defined $file2;
    my $id_seq = read_fasta($fasta);
    my @ids = sort keys %$id_seq;
    die "Error: FASTA file must contain exactly 2 sequences (found "
        . scalar(@ids) . ").\n" unless @ids == 2;
    ($id1, $id2) = @ids;
    $seq1 = uc $id_seq->{$id1};
    $seq2 = uc $id_seq->{$id2};
}
elsif (defined $seq1_in && defined $seq2_in) {
    die "Cannot mix -s1/-s2 with -f/-f1/-f2.\n"
        if defined $fasta || defined $file1 || defined $file2;
    ($id1, $id2) = ('seq1', 'seq2');
    $seq1 = uc $seq1_in;
    $seq2 = uc $seq2_in;
}
elsif (defined $file1 && defined $file2) {
    die "Cannot mix -f1/-f2 with -f/-s1/-s2.\n"
        if defined $fasta || defined $seq1_in || defined $seq2_in;
    ($id1, $id2) = ('seq1', 'seq2');
    $seq1 = uc read_single_line($file1);
    $seq2 = uc read_single_line($file2);
}
else {
    die $usage;
}

# Validate characters (allow letters only; case-insensitive, uppercased above)
for my $s ($seq1, $seq2) {
    die "Error: sequence contains non-letter characters: $s\n"
        unless $s =~ /\A[A-Za-z]*\z/;
}

# Auto-detect sequence type
if ($type eq 'auto') {
    my $combined = $seq1 . $seq2;
    if ($combined =~ /[^ACGTUN]/i) {
        $type = 'protein';
    }
    else {
        $type = 'dna';
    }
}

# Resolve scoring scheme
my $sub_matrix;    # hashref for substitution matrix (undef = simple mode)
if ($type eq 'protein' && $matrix eq 'blosum62') {
    $sub_matrix = get_blosum62();
    $gap = -8 unless defined $gap;
}
else {
    # Simple match/mismatch mode
    if ($type eq 'protein') {
        $match   =  5 unless defined $match;
        $unmatch = -4 unless defined $unmatch;
        $gap     = -8 unless defined $gap;
    }
    else {
        $match   =  1 unless defined $match;
        $unmatch = -1 unless defined $unmatch;
        $gap     = -2 unless defined $gap;
    }
}

# Redirect output if requested
if (defined $output) {
    open my $out, '>', $output or die "Cannot write to $output: $!\n";
    select $out;
}

#==============================================================================
# Run Needleman-Wunsch alignment
#==============================================================================
my $t0 = join('.', gettimeofday);
my ($score, $aln1, $aln2, $path, $score_matrix) = needleman_wunsch($seq1, $seq2,
                                                    $match, $unmatch, $gap,
                                                    $sub_matrix);
my $t1 = join('.', gettimeofday);

# Build match line and identity ratio
my $match_line = '';
for my $k (0 .. length($aln1) - 1) {
    $match_line .= substr($aln1, $k, 1) eq substr($aln2, $k, 1) ? '|' : ' ';
}
my $match_count = ($match_line =~ tr/|/|/);
my $identity    = length($aln1) ? ($match_count / length($aln1)) * 100 : 0;

# Verbose: dump score matrix and path (only in detailed mode)
if ($verbose && $detailed) {
    print "----------------------------------------------------------------------\n";
    print "Type: $type | Matrix: $matrix | Gap: $gap\n";
    if ($sub_matrix) {
        print "Scoring: BLOSUM62 substitution matrix\n";
    }
    else {
        print "Scoring: match=$match mismatch=$unmatch gap=$gap\n";
    }
    print "----------------------------------------------------------------------\n";
    print "Score matrix (rows = seq1, cols = seq2):\n";
    my @s1 = split //, $seq1;
    my @s2 = split //, $seq2;
    printf "%6s %6s", '', '-';
    for my $c (@s2) { printf "%6s", $c }
    print "\n";
    for my $i (0 .. scalar(@s1)) {
        printf "%6s", $i == 0 ? '-' : $s1[$i - 1];
        for my $j (0 .. scalar(@s2)) {
            printf "%6d", $score_matrix->[$i][$j];
        }
        print "\n";
    }
    print "\nBacktrack path (H=gap in seq2, V=gap in seq1, D=diagonal):\n";
    print $path, "\n\n";
}

# Output modes
if ($detailed) {
    print "----------------------------------------------------------------------\n";
    printf "Sequence type      : %s\n", $type;
    printf "Scoring scheme     : %s\n", ($sub_matrix ? "BLOSUM62" : "simple ($match/$unmatch/$gap)");
    printf "Sequence 1 (%s), length: %d\n", $id1, length($seq1);
    printf "Sequence 2 (%s), length: %d\n", $id2, length($seq2);
    print  "----------------------------------------------------------------------\n";
    printf "Alignment score    : %d\n", $score;
    printf "Identity ratio     : %.2f%%  (%d/%d)\n", $identity, $match_count, length($aln1);
    printf "Aligned length     : %d\n", length($aln1);
    print  "----------------------------------------------------------------------\n";
    print  "Aligned seq1 : $aln1\n";
    print  "               $match_line\n";
    print  "Aligned seq2 : $aln2\n";
    print  "----------------------------------------------------------------------\n";
    printf "Run time           : %.4fs\n", $t1 - $t0;
    print  "----------------------------------------------------------------------\n";
}

# FASTA output (always output, regardless of -d)
print ">$id1\n$aln1\n";
print ">$id2\n$aln2\n";

#==============================================================================
# Subroutines
#==============================================================================

# Read a FASTA file; returns hashref { id => concatenated_sequence }
sub read_fasta {
    my $file = shift;
    open my $in, '<', $file or die "Cannot open $file: $!\n";
    my (%seq, $id);
    while (<$in>) {
        chomp;
        s/\r$//;
        next unless length;
        if (/^>(.*)/) {
            $id = (split(/\s+/, $1))[0];
        }
        elsif (defined $id) {
            $seq{$id} .= uc $_;
        }
    }
    close $in;
    die "Error: no FASTA records found in $file.\n" unless keys %seq;
    return \%seq;
}

# Read first non-empty line from a plain text file (sequence only, no header)
sub read_single_line {
    my $file = shift;
    open my $in, '<', $file or die "Cannot open $file: $!\n";
    while (<$in>) {
        chomp;
        s/\r$//;
        next unless length;
        close $in;
        return $_;
    }
    close $in;
    die "Error: no sequence found in $file.\n";
}

# Core Needleman-Wunsch algorithm
#   $sub_matrix: hashref { 'AA' => 4, 'AR' => -1, ... } or undef for simple mode
# Returns: ($best_score, $aligned_seq1, $aligned_seq2, $path_string, \@score_matrix)
# Path chars: 'D' (diagonal), 'H' (gap in seq2, horizontal), 'V' (gap in seq1, vertical)
sub needleman_wunsch {
    my ($s1, $s2, $m, $u, $g, $sub_matrix) = @_;
    my @a1 = split //, $s1;
    my @a2 = split //, $s2;
    my $n1 = scalar @a1;
    my $n2 = scalar @a2;

    # Edge case: empty sequence(s)
    if ($n1 == 0 && $n2 == 0) {
        return (0, '', '', '', []);
    }
    if ($n1 == 0) {
        my @sc;
        $sc[0][$_] = $_ * $g for 0 .. $n2;
        return ($n2 * $g, '-' x $n2, $s2, 'V' x $n2, \@sc);
    }
    if ($n2 == 0) {
        my @sc;
        $sc[$_][0] = $_ * $g for 0 .. $n1;
        return ($n1 * $g, $s1, '-' x $n1, 'H' x $n1, \@sc);
    }

    # Initialize score and arrow matrices
    # arrow encoding: 0 = diagonal, 1 = up, 2 = left
    my @score;
    my @arrow;
    $score[0][$_] = $_ * $g for 0 .. $n2;
    $arrow[0][$_] = 2     for 0 .. $n2;
    for my $i (1 .. $n1) {
        $score[$i][0] = $i * $g;
        $arrow[$i][0] = 1;
    }

    # Fill matrices
    for my $i (1 .. $n1) {
        for my $j (1 .. $n2) {
            my $diag;
            if ($sub_matrix) {
                my $key = $a1[$i-1] . $a2[$j-1];
                $diag = $score[$i-1][$j-1] + (exists $sub_matrix->{$key} ? $sub_matrix->{$key} : -4);
            }
            else {
                $diag = $score[$i-1][$j-1] + ($a1[$i-1] eq $a2[$j-1] ? $m : $u);
            }
            my $up   = $score[$i-1][$j]  + $g;
            my $left = $score[$i][$j-1]  + $g;

            # Tie priority: diagonal > up > left
            if ($diag >= $up && $diag >= $left) {
                $score[$i][$j] = $diag;
                $arrow[$i][$j] = 0;
            }
            elsif ($up >= $left) {
                $score[$i][$j] = $up;
                $arrow[$i][$j] = 1;
            }
            else {
                $score[$i][$j] = $left;
                $arrow[$i][$j] = 2;
            }
        }
    }

    # Traceback from (n1, n2) to (0, 0)
    my ($i, $j) = ($n1, $n2);
    my ($aln1, $aln2, $path) = ('', '', '');
    while ($i > 0 || $j > 0) {
        my $a = $arrow[$i][$j];
        if ($a == 0) {           # diagonal
            $aln1 .= $a1[$i-1];
            $aln2 .= $a2[$j-1];
            $path .= 'D';
            $i--; $j--;
        }
        elsif ($a == 1) {        # up (gap in seq2)
            $aln1 .= $a1[$i-1];
            $aln2 .= '-';
            $path .= 'H';
            $i--;
        }
        else {                   # left (gap in seq1)
            $aln1 .= '-';
            $aln2 .= $a2[$j-1];
            $path .= 'V';
            $j--;
        }
    }
    $aln1 = reverse $aln1;
    $aln2 = reverse $aln2;
    $path = reverse $path;

    return ($score[$n1][$n2], $aln1, $aln2, $path, \@score);
}

# BLOSUM62 substitution matrix (24x24: 20 standard AA + B, Z, X, *)
# Returns hashref { 'AR' => -1, 'AA' => 4, ... }
sub get_blosum62 {
    my @header = qw(A R N D C Q E G H I L K M F P S T W Y V B Z X *);
    my @rows = (
        # A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
        [  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -4 ], # A
        [ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4 ], # R
        [ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4 ], # N
        [ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4 ], # D
        [  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4 ], # C
        [ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4 ], # Q
        [ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4 ], # E
        [  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4 ], # G
        [ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4 ], # H
        [ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4 ], # I
        [ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4 ], # L
        [ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4 ], # K
        [ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4 ], # M
        [ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4 ], # F
        [ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -1, -4 ], # P
        [  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4 ], # S
        [  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4 ], # T
        [ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4 ], # W
        [ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4 ], # Y
        [  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4 ], # V
        [ -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4 ], # B
        [ -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4 ], # Z
        [ -1, -1, -1, -1, -2, -1, -1, -2, -1, -1, -1, -1, -1, -2, -1,  0,  0, -2, -1, -1, -1, -1, -1, -4 ], # X
        [ -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1 ], # *
    );

    my %matrix;
    for my $i (0 .. $#header) {
        for my $j (0 .. $#header) {
            $matrix{$header[$i] . $header[$j]} = $rows[$i][$j];
        }
    }
    return \%matrix;
}
