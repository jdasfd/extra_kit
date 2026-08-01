#!/usr/bin/perl -w
#
# fasta2axt.pl -- Convert FASTA alignment (2 sequences) to AXT format
#                 (KaKs_Calculator compliant)
#
# Author: Yuqian Jiang
# Created: 2026-08-02
# Supported: GLM 5
# Referenced:
#   - parseFastaIntoAXT.pl (KaKs_Calculator, Zhang Z et al. 2006)
# Change logs
# 2026-08-02: refactor as parameterized converter; modular design; fix AXT format to match KaKs_Calculator specification
#

use strict;
use warnings;
use Getopt::Long;

my ($fasta, $output, $help);

GetOptions(
    'fasta|f=s'  => \$fasta,
    'output|o=s' => \$output,
    'help|h+'    => \$help,
) or die "Use -h for help.\n";

my $usage = <<'__GUIDE__';
Convert FASTA alignment (exactly 2 sequences) to AXT format

Usage:
    perl fasta2axt.pl -f <input.fasta>                       [options]
    perl fasta2axt.pl -f <input.fasta> -o <output.axt>
    cat aln.fasta | perl fasta2axt.pl -f -

Required:
    -f|--fasta  <file> : input FASTA file (or '-' for stdin)

Optional:
    -o|--output <file> : output AXT file (default: STDOUT)
    -h|--help          : show this help

__GUIDE__

die $usage if $help;
die $usage unless defined $fasta;

# Redirect output if requested
if (defined $output && $output ne '-') {
    open my $out, '>', $output or die "Cannot write to $output: $!\n";
    select $out;
}

#==============================================================================
# Parse FASTA
#==============================================================================
my $records = parse_fasta($fasta);

# Validate: exactly 2 sequences
die "Error: Input must contain exactly 2 sequences for AXT format. Found "
    . scalar(@$records) . "\n" unless @$records == 2;

my ($id1, $seq1) = @{$records->[0]};
my ($id2, $seq2) = @{$records->[1]};

# Validate equal length
if (length($seq1) != length($seq2)) {
    warn "Warning: Sequences are not equal in length! "
       . "Seq1: " . length($seq1) . ", Seq2: " . length($seq2) . "\n"
       . "KaKs_Calculator will likely fail. Please check your alignment.\n";
}

#==============================================================================
# Write AXT (KaKs_Calculator format)
#==============================================================================
my $caption = "$id1-$id2";
print "$caption\n";
print "$seq1\n";
print "$seq2\n";
print "\n";    # AXT blocks are separated by a blank line

#==============================================================================
# Subroutines
#==============================================================================

# Parse a FASTA file (or stdin if $file eq '-')
# Returns: arrayref of [id, sequence] pairs (sequence uppercased, whitespace stripped)
#   - id is the first whitespace-separated token after '>'
sub parse_fasta {
    my $file = shift;
    my $fh;
    if ($file eq '-') {
        $fh = \*STDIN;
    }
    else {
        open $fh, '<', $file or die "Cannot open $file: $!\n";
    }

    my @records;
    my ($id, $seq);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        if ($line =~ /^>(\S+)/) {
            if (defined $id) {
                push @records, [$id, $seq];
            }
            $id  = $1;
            $seq = '';
        }
        elsif (defined $id) {
            $line =~ s/\s+//g;
            $seq .= uc $line;
        }
    }
    if (defined $id) {
        push @records, [$id, $seq];
    }
    close $fh unless $file eq '-';
    die "Error: no FASTA records found in $file.\n" unless @records;
    return \@records;
}
