#!/usr/bin/perl -w
#
# domain_arrange.pl -- Filtering and rearranging domains inside a protein
#
# Author: Yuqian Jiang
# Created: 2025-04-22

use strict;
use warnings;
use AlignDB::IntSpan;
use GetOpt::Long;

=head1 NAME

domain_arrange.pl - Filtering and rearranging domains inside a protein

=head1 SYNOPSIS

    domain_arrange.pl
        Filtering and rearranging domains inside a protein

    Usage:
        perl domain_arrange.pl -i <domain_scan.tsv>

    Note:
        .tsv need satisfied the following rules for usage
        Col1: Gene/Protein ids
        Col2: Domain
        Col3: Domain start
        Col4: Domain end
        Col5: E-value
        Print to stdout, please redirect to the file

=cut

GetOptions(
    'input|i=s'     => \(my $input),
    "h|help"        => sub { Getopt::Long::HelpMessage(0) }
) or Getopt::Long::HelpMessage(1);

if ( !defined $input ) {
    print STDERR "Error: cannot find input tsv.\n";
    die Getopt::Long::HelpMessage(1);
}

#----------------------------------------------------------#
# Main program
#----------------------------------------------------------#

my $gene;
my @print;
my $domain_set = AlignDB::IntSpan -> new;

open my $TSV_IN, "<", $input or die "Error: cannot read file in.\n";

while (<$TSV_IN>) {
    chomp;
    my $print_list = $_;
    my @array = split/\t/, $_;
    if ( ! defined $gene ) {
        $gene = $array[0];
        push @print, $print_list;
        $domain_set -> add_range($array[2], $array[3]);
    }
    else {
        if ( $array[0] eq $gene ) {
            my $test_set = AlignDB::IntSpan -> new;
            $test_set -> add_range($array[2], $array[3]);
            my $result = $domain_set -> intersect($test_set);
            if ( $result -> is_empty  ){
                $domain_set -> add_range($array[2], $array[3]);
                push @print, $print_list;
            }
        }
        else {
            print join "\n", @print;
            print "\n";
            $domain_set -> clear;
            $gene = $array[0];
            @print = ();
            push @print, $print_list;
            $domain_set -> add_range($array[2], $array[3])
        }
    }
}

END {
    print join "\n", @print;
    print "\n";
}
