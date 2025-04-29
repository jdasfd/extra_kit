#!/usr/bin/perl -w
#
# domain_arrange.pl -- Filtering and rearranging domains inside a protein
#
# Author: Yuqian Jiang
# Created: 2023-03-22
# Updated: 2025-04-29 Add more sub programs

use strict;
use warnings;
use AlignDB::IntSpan;
use Getopt::Long;
use Math::BigFloat;

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

our %DOMAIN;
Math::BigFloat -> accuracy(50);

load_domains($input, \%DOMAIN);

for my $gene (sort keys %DOMAIN){
    my @sorted = sort_domains($DOMAIN{$gene});
    my @filtered = filter_domains(@sorted);
    @filtered = sort { $a -> {start} <=> $b -> {start} } @filtered;
    for (@filtered) {
        print $_ -> {line} . "\n";
    }
}

#----------------------------------------------------------#
# Sub program
#----------------------------------------------------------#

sub load_domains {
    my ($file, $hash) = @_;
    open my $fh, "<", $file or die "Cannot open $file: $!";

    while (<$fh>) {
        chomp;
        my @fields = split /\t/;
        next unless @fields == 5;

        my ($gene, $name, $start, $end, $evalue) = @fields;
        $evalue = Math::BigFloat -> new($evalue);

        my $index = scalar @{ $hash -> {$gene} // [] };
        push @{ $hash->{$gene} }, {
            index  => $index,
            name   => $name,
            start  => $start,
            end    => $end,
            evalue => $evalue,
            line   => $_
        };
    }
    close $fh;
}

sub sort_domains {
    my ($domains) = @_;
    return sort {
        $a->{evalue} <=> $b->{evalue} ||
        $a->{index} <=> $b->{index}
    } @$domains;
}

sub filter_domains {
    my @domains = @_;
    my @result;
    my $span = AlignDB::IntSpan -> new;

    foreach my $dom (@domains) {
        my $test = AlignDB::IntSpan -> new;
        $test -> add_range($dom -> {start}, $dom -> {end});

        if ($span -> intersect($test) -> is_empty) {
            $span -> add($test);
            push @result, $dom;
        }
    }
    return @result;
}

__END__
