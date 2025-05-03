#!/usr/bin/perl -w
#
# extract_longest_ncbi.pl -- extract longest peptides or proteins from ncbi gff3
#
# Author: Yuqian Jiang
# Created: 2025-05-03

use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;
use Data::Dumper;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

extract_longest_ncbi.pl - extract longest peptides or proteins from ncbi gff3

=head1 SYNOPSIS

    perl extract_longest_ncbi.pl -i <gff3> -o <output>
      Options:
        --gff           -g  STR     gff annotation file
        --output        -o  STR     output files (also as the gff3 format), default: STDOUT
        --help          -h          brief help message

=cut

GetOptions(
    'help|h'     => sub { Getopt::Long::HelpMessage(0) },
    'gff|g=s'    => \( my $gff_file ),
    'output|o=s' => \( my $output = 'stdout' )
) or Getopt::Long::HelpMessage(1);

if ( !defined $gff_file ) {
    print STDERR "Error: please supply an annotation gff3 format file.\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($gff_file) -> is_file ) {
    die "Error: can't find file [$gff_file].\n";
}

#----------------------------------------------------------#
# Main program
#----------------------------------------------------------#

my (%MRNA_GENE, %GENE_FORMAT);
my $line_num = 0;

# read gff3 in
open my $GFF_IN, "<", $gff_file;

while (<$GFF_IN>) {
    $line_num++;
    chomp;
    next if (/^\#/ || /^\s+$/ || /^$/);

    my @arrays = split/\t/;
    if ($arrays[2] eq "gene") {
        my ($gene_id) = $arrays[8] =~ /ID=([^;]+)/;
        $GENE_FORMAT{$gene_id} = {};
    }
    elsif ($arrays[2] eq "mRNA" || $arrays[2] eq "transcript") {
        my ($transcript_id) = $arrays[8] =~ /ID=([^;]+)/;
        my ($gene_parent) = $arrays[8] =~ /Parent=([^;]+)/;
        unless (exists $GENE_FORMAT{$gene_parent}) {
            die "Error at line $line_num: can't find gene parent, index error!\n";
        }
        else {
            $MRNA_GENE{$transcript_id} = $gene_parent;
            $GENE_FORMAT{$gene_parent}{$transcript_id} = 0;
        }
    }
    elsif ($arrays[2] eq "CDS") {
        my ($transcript_parent) = $arrays[8] =~ /Parent=([^;]+)/;
        my $length = $arrays[4] - $arrays[3] + 1;
        if (my $gene_id = $MRNA_GENE{$transcript_parent}) {
            $GENE_FORMAT{$gene_id}{$transcript_parent} += $length;
        }
    }
}

# output
my $out_fh;
if ( lc($output) eq "stdout" ) {
    $out_fh = *STDOUT;
}
else {
    open $out_fh, ">", $output;
}

foreach my $gene_id (keys %GENE_FORMAT) {
    my ($max_len, $longest_mrna) = (0, '');
    while (my ($mrna, $len) = each %{$GENE_FORMAT{$gene_id}}) {
        ($max_len, $longest_mrna) = ($len, $mrna) if $len > $max_len;
    }
    next if $longest_mrna eq ""; # skip those ncrna
    print "$gene_id\t$longest_mrna\t$max_len\n";
}

close $out_fh;
