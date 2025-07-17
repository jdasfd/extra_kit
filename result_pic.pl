#!/usr/bin/perl -w
#
# A simple script for generating a matrix.
#
# Author: Yuqian Jiang
# Created: 2023-03-22
# Decided to put it into RAID tool-box for analyzing

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Path::Tiny;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

result_pic.pl - a simple script for show elements existence or their repeated times

=head1 SYNOPSIS

    perl result_pic.pl -l <list_file> -n <headline_name> -i <tsv_file> -t <circle|count>
      Options:
        --help          -h                  brief help message
        --list          -l  STR             list file contains elements for counting
        --name          -n  STR             headline first name in (A\B) format
        --input         -i  STR             counting files in tsv (col1: sample names, col2: elements)
        --type          -t  circle|count    output circle (contains or not) or count (contains how many)

=cut

GetOptions(
    'help|h'        => sub { Getopt::Long::HelpMessage(0) },
    'list|l=s'      => \( my $list_file ),
    'name|n=s'      => \( my $name ),
    'input|i=s'     => \( my $in_file ),
    'type|t=s'      => \( my $type ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $list_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($list_file)->is_file ) {
    die "Error: can't find file [$list_file]";
}

if ( !defined $in_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($in_file)->is_file ) {
    die "Error: can't find file [$in_file]";
}

if ( !defined $type && $type ne "circle" && $type ne "count" ) {
    die "Please specify your output type: circle/count";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

# variable
our (%COUNT, %SAMPLE);
our (@headline, @order);
my $element = "";

# dealing with list input
open my $LI_IN, '<', $list_file;
while (<$LI_IN>) {
    chomp;
    $COUNT{$_} = 0;
    push @headline, $_;
}

close $LI_IN;

open my $DEAL_IN, '<', $in_file;
while (<$DEAL_IN>) {
    chomp;
    my @array = split/\t/, $_;
    if ($element ne "$array[0]") {
        push @order, $array[0];
        $element = $array[0];
    }
    push @{$SAMPLE{$array[0]}}, $array[1];
}

close $DEAL_IN;

# create headline
&create_head($name, @headline);

# count for every sample
foreach my $order ( @order ) {
    for my $count ( @{$SAMPLE{$order}} ) {
        $COUNT{$count}++;
    }
    print "$order\t";
    my @printlist;
    if ( $type eq "circle" ) {
        foreach my $head ( @headline ) {
            if ( $COUNT{$head} != 0 ) {
                $COUNT{$head} = "o";
            }
            else {
                $COUNT{$head} = "x";
            }
            push @printlist, $COUNT{$head};
        }
    }
    elsif ( $type eq "count" ) {
        foreach my $head ( @headline ) {
            push @printlist, $COUNT{$head};
        }
    }
    print join ("\t", @printlist);
    print "\n";
    # empty hash COUNT
    for my $key ( keys %COUNT ) {
        $COUNT{$key} = 0;
    }
}

#----------------------------------------------------------#
# sub to create headline
#----------------------------------------------------------#

sub create_head {
    my ($scalar, @list) = @_;
    print join ("\t", @_);
    print "\n";
}

__END__
