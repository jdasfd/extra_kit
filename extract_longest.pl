#!/usr/bin/perl -w
#
# extract_longest.pl -- extract longest CDS transcripts per gene from gff3
#
# Author: Yuqian Jiang
# Created: 2025-05-02
# Updated: 2025-05-03 Deal with the pseudogene in the script
# Updated: 2025-05-03 Filter the one without CDS regions
# Updated: 2026-09-01 Add hash sort feature (by transcript id) to avoid random output
# Updated: 2026-09-03 Rename from extract_longest_ncbi.pl and add --ncbi|--liftoff option
# Updated: 2026-09-03 Liftoff diagnostics go to a separate log (--log, default STDERR), table stays 4 columns
# Updated: 2026-09-03 Add [single_transcript] gene-level log entry; one gene-level tag per gene for direct counting

use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;
use Data::Dumper;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

extract_longest.pl - extract longest CDS transcripts per gene from gff3

=head1 SYNOPSIS

    perl extract_longest.pl (--ncbi|--liftoff) -i <gff3> -o <output>
      Options:
        --ncbi              gff3 comes from NCBI (original logic, untouched)
        --liftoff           gff3 comes from the Liftoff/AGAT pipeline; isoforms with
                            CDS always win over those without (among CDS ones the
                            longest CDS is kept, AGAT-like but order-independent);
                            every gene gets exactly one gene-level log entry, plus
                            isoform-level [no_cds] entries:
                                [single_transcript]  gene with only one mRNA
                                [multi_mrna]         gene with multiple mRNAs (with CDS lengths)
                                [no_cds]             isoforms without CDS
                                [orphan]             mRNA whose parent gene is not found
        --gff           -g  STR     gff annotation file
        --output        -o  STR     output files (also as the gff3 format), default: STDOUT
        --log           -l  STR     diagnostic log file (liftoff mode only), default: STDERR
        --help          -h          brief help message

=cut

GetOptions(
    'help|h'     => sub { Getopt::Long::HelpMessage(0) },
    'ncbi'       => \( my $flag_ncbi ),
    'liftoff'    => \( my $flag_liftoff ),
    'gff|g=s'    => \( my $gff_file ),
    'output|o=s' => \( my $output = 'stdout' ),
    'log|l=s'    => \( my $log_file )
) or Getopt::Long::HelpMessage(1);

my $source_flag_num = ( defined $flag_ncbi ? 1 : 0 ) + ( defined $flag_liftoff ? 1 : 0 );
if ( $source_flag_num != 1 ) {
    print STDERR "Error: please declare the annotation source by exactly one of --ncbi or --liftoff.\n";
    die Getopt::Long::HelpMessage(1);
}
my $is_liftoff = defined $flag_liftoff ? 1 : 0;

if ( !defined $gff_file ) {
    print STDERR "Error: please supply an annotation gff3 format file.\n";
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($gff_file) -> is_file ) {
    die "Error: can't find file [$gff_file].\n";
}

# log filehandle (liftoff diagnostics)
my $log_fh = *STDERR;
if ( defined $log_file ) {
    if ($is_liftoff) {
        open $log_fh, ">", $log_file or die "Can't write to $log_file: $!";
    }
    else {
        print STDERR "Warning: --log only applies to --liftoff mode, option ignored.\n";
    }
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
    if ($arrays[2] =~ /^(gene|pseudogene)$/) {
        my ($gene_id) = $arrays[8] =~ /ID=([^;]+)/;
        $GENE_FORMAT{$gene_id} = { type => $arrays[2], mrnas => {}};
    }
    elsif ($arrays[2] =~ /^(mRNA|transcript)$/) {
        my ($transcript_id) = $arrays[8] =~ /ID=([^;]+)/;
        my ($gene_parent) = $arrays[8] =~ /Parent=([^;]+)/;

        # liftoff: log the mRNA whose parent gene is not indexed
        if ( $is_liftoff && !exists $GENE_FORMAT{$gene_parent} ) {
            print $log_fh "[orphan]\tline $line_num\tmrna="
                . ( defined $transcript_id ? $transcript_id : '--' )
                . "\tparent=" . ( defined $gene_parent ? $gene_parent : '--' ) . "\n";
            next;
        }
        unless (exists $GENE_FORMAT{$gene_parent}) {
            die "Error at line $line_num: Paret $gene_parent not found, index error!\n";
        }

        $MRNA_GENE{$transcript_id} = $gene_parent;
        $GENE_FORMAT{$gene_parent} -> {mrnas} -> {$transcript_id} = 0;
    }
    elsif ($arrays[2] eq "CDS") {
        my ($transcript_parent) = $arrays[8] =~ /Parent=([^;]+)/;
        next unless exists $MRNA_GENE{$transcript_parent};

        my $length = $arrays[4] - $arrays[3] + 1;
        my $gene_id = $MRNA_GENE{$transcript_parent};
        $GENE_FORMAT{$gene_id} -> {mrnas} {$transcript_parent} += $length;
    }
}


# output
my $out_fh;
if ( lc($output) eq "stdout" ) {
    $out_fh = *STDOUT;
}
else {
    open $out_fh, ">", $output or die "Can't write to $output: $!";
}

print $out_fh "Gene_ID\tLongest_mRNA\tLength\tType\n";

for my $gene_id (sort keys %GENE_FORMAT) {
    my $gene_info = $GENE_FORMAT{$gene_id};
    my $type = $gene_info -> {type};
    my %mrnas = %{$gene_info -> {mrnas}};

    next unless %mrnas;

    # liftoff: audit the isoforms of this gene in the log
    if ($is_liftoff) {
        my @mrna_ids = sort keys %mrnas;
        my $gene_tag = @mrna_ids > 1 ? "multi_mrna" : "single_transcript";
        print $log_fh "[$gene_tag]\t$gene_id\t" . scalar(@mrna_ids) . "\t"
            . join(";", map {"$_=$mrnas{$_}"} @mrna_ids) . "\n";
        for my $mrna_id (@mrna_ids) {
            print $log_fh "[no_cds]\t$gene_id\t$mrna_id\n" if $mrnas{$mrna_id} == 0;
        }
    }

    my ($longest, $max_len) = ('', 0);
    # isoforms with CDS always win over those without (CDS-less length is 0);
    # among CDS-containing ones keep the longest CDS, tie -> smallest mRNA ID
    for my $mrna (sort keys %mrnas) {
        my $len = $mrnas{$mrna};
        ($longest, $max_len) = ($mrna, $len) if $len > $max_len;
    }
    next if $max_len == 0;

    print $out_fh "$gene_id\t$longest\t$max_len\t$type\n";
}

close $out_fh;
close $log_fh if fileno($log_fh) > 2;
