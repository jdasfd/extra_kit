#!/usr/bin/perl -w
#
# calculate_4DTv.pl -- Calculate 4DTv (transversion rate on 4-fold degenerate sites) from codon alignments in AXT format
#
# Author: Yuqian Jiang (refactored from original)
# Original authors: Sun Ming'an, Fan Wei, Li Jun (BGI, 2008)
# Change logs
# 2026-07-31: refactor as parameterized calculator; add multiple codon table support; auto-detect 4-fold degenerate codons; modular design

use strict;
use warnings;
use Getopt::Long;

my ($axt, $output, $table, $help);

GetOptions(
    'axt|a=s'    => \$axt,
    'output|o=s' => \$output,
    'table|t=s'  => \$table,
    'help|h+'    => \$help,
) or die "Use -h for help.\n";

my $usage = <<'__GUIDE__';
    Calculate 4DTv (transversion rate at 4-fold degenerate sites) from codon
    alignments in AXT format.

Usage:
    perl calculate_4DTv.pl -a <input.axt>                       [options]
    perl calculate_4DTv.pl -a <input.axt> -o <output.tsv>
    cat in.axt | perl calculate_4DTv.pl -a -

Required:
    -a|--axt   <file> : input AXT file (or '-' for stdin)

Optional:
    -o|--output<file> : output TSV file (default: STDOUT)
    -t|--table  <str> : codon table:
                        standard     (NCBI 1,  default)
                        vert-mito    (NCBI 2)
                        inverte-mito (NCBI 5)
                        plant-plastid(NCBI 11)
    -h|--help         : show this help

Reference:
    Hasegawa M, Kishino H, Yano T. J. Mol. Evol. 22(2):160, 1985.
__GUIDE__

die $usage if $help;
die $usage unless defined $axt;

$table = 'standard' unless defined $table;
$table = lc $table;

# Resolve codon table and 4-fold degenerate codons
my $codon_table   = get_codon_table($table);
my %codons        = %$codon_table;
my %four_fold     = get_four_fold_codons($codon_table);

# Redirect output if requested
if (defined $output && $output ne '-') {
    open my $out, '>', $output or die "Cannot write to $output: $!\n";
    select $out;
}

#==============================================================================
# Parse AXT and compute 4DTv for each block
#==============================================================================
my $blocks = parse_axt($axt);

print "tag\t4dtv_corrected\t4dtv_raw\tcodon_4d\tcodon_4dt\n";
for my $block (@$blocks) {
    my ($tag, $seq1, $seq2) = @$block;
    my ($corrected_4dtv, $raw_4dtv, $codon_4d, $codon_4dt) =
        calculate_4dtv($seq1, $seq2, \%codons, \%four_fold);
    print "$tag\t$corrected_4dtv\t$raw_4dtv\t$codon_4d\t$codon_4dt\n";
}

#==============================================================================
# Subroutines
#==============================================================================

# Parse AXT file (or stdin)
# Returns: arrayref of [tag, seq1, seq2]
sub parse_axt {
    my $file = shift;
    my $fh;
    if ($file eq '-') {
        $fh = \*STDIN;
    }
    else {
        open $fh, '<', $file or die "Cannot open $file: $!\n";
    }

    local $/ = "\n\n";    # block separator
    my @blocks;
    while (my $block = <$fh>) {
        chomp $block;
        $block =~ s/\r$//;
        next unless $block =~ /\S/;
        if ($block =~ /^(\S+)\n(\S+)\n(\S+)/) {
            push @blocks, [$1, $2, $3];
        }
        else {
            warn "Warning: malformed AXT block skipped.\n";
        }
    }
    close $fh unless $file eq '-';
    die "Error: no AXT blocks found in $file.\n" unless @blocks;
    return \@blocks;
}

# Calculate 4DTv using HKY85 model
#   $codons_ref     : hashref codon → amino acid
#   $four_fold_ref  : hashref of 4-fold degenerate codons (codon => 1)
# Returns: ($corrected_4dtv, $raw_4dtv, $codon_4d, $codon_4dt)
sub calculate_4dtv {
    my ($str1, $str2, $codons_ref, $four_fold_ref) = @_;
    my %codons     = %$codons_ref;
    my %four_fold  = %$four_fold_ref;

    my ($codon_4d, $codon_4dt) = (0, 0);
    my %fre = ();    # base frequency at 4-fold degenerate sites
    my %transversion = (
        A => 'TC', C => 'AG', G => 'TC', T => 'AG',
    );

    for (my $i = 0; $i < length($str1); $i += 3) {
        my $codon1 = uc substr($str1, $i, 3);
        my $codon2 = uc substr($str2, $i, 3);

        # Both codons must be 4-fold degenerate and encode the same amino acid
        next unless exists $four_fold{$codon1} && exists $four_fold{$codon2};
        next unless $codons{$codon1} eq $codons{$codon2};

        my $base1 = substr($codon1, 2, 1);
        my $base2 = substr($codon2, 2, 1);

        $fre{$base1}++;
        $fre{$base2}++;
        $codon_4d++;
        $codon_4dt++ if is_transversion($base1, $base2, \%transversion);
    }

    my ($V, $d);
    if ($codon_4d > 0) {
        $V = $codon_4dt / $codon_4d;    # raw 4DTv

        # HKY85 correction
        $fre{Y} = ($fre{T} // 0) + ($fre{C} // 0);
        $fre{R} = ($fre{A} // 0) + ($fre{G} // 0);
        for my $k (keys %fre) {
            $fre{$k} = 0.5 * $fre{$k} / $codon_4d;
        }

        if (   $fre{Y} != 0 && $fre{R} != 0
            && $fre{A} != 0 && $fre{C} != 0
            && $fre{G} != 0 && $fre{T} != 0) {

            my $a = -1 * log(
                1 - $V * ($fre{T}*$fre{C}*$fre{R}/$fre{Y}
                         + $fre{A}*$fre{G}*$fre{Y}/$fre{R})
                / (2 * ($fre{T}*$fre{C}*$fre{R} + $fre{A}*$fre{G}*$fre{Y}))
            );

            if (1 - $V / (2 * $fre{Y} * $fre{R}) > 0) {
                my $b = -1 * log(1 - $V / (2 * $fre{Y} * $fre{R}));
                $d = 2 * $a * ($fre{T}*$fre{C}/$fre{Y} + $fre{A}*$fre{G}/$fre{R})
                   - 2 * $b * ($fre{T}*$fre{C}*$fre{R}/$fre{Y}
                              + $fre{A}*$fre{G}*$fre{Y}/$fre{R}
                              - $fre{Y}*$fre{R});
            }
            else {
                $d = 'NA';
            }
        }
        else {
            $d = 'NA';
        }
    }
    else {
        $V = 'NA';
        $d = 'NA';
    }

    return ($d, $V, $codon_4d, $codon_4dt);
}

# Check if base1 → base2 is a transversion
sub is_transversion {
    my ($base1, $base2, $transversion_ref) = @_;
    return 0 unless exists $transversion_ref->{$base1};
    return $transversion_ref->{$base1} =~ /$base2/ ? 1 : 0;
}

# Auto-detect 4-fold degenerate codons from a codon table
# A codon family is "4-fold degenerate" if all four codons sharing the same
# first two nucleotides encode the same amino acid.
# Returns: hashref (codon => 1)
sub get_four_fold_codons {
    my $codon_table = shift;
    my @bases = qw(A C G T);
    my %four_fold;

    # Group codons by first two nucleotides (XY*)
    my %prefix_aa;
    for my $b1 (@bases) {
        for my $b2 (@bases) {
            my $prefix = $b1 . $b2;
            my @aas;
            my @codons;
            for my $b3 (@bases) {
                my $codon = $prefix . $b3;
                if (exists $codon_table->{$codon}) {
                    push @aas,    $codon_table->{$codon};
                    push @codons, $codon;
                }
                else {
                    # incomplete table - skip this family
                    @aas = ();
                    last;
                }
            }
            # All four codons must encode the same amino acid
            if (@aas == 4) {
                my $first = $aas[0];
                my $all_same = 1;
                for my $aa (@aas) {
                    if ($aa ne $first) { $all_same = 0; last; }
                }
                if ($all_same) {
                    $four_fold{$_} = 1 for @codons;
                }
            }
        }
    }
    return %four_fold;
}

# Return codon table as hashref
# Each sub returns hashref codon → amino acid (single letter)
sub get_codon_table {
    my $name = shift;

    my %tables = (
        standard      => \&_table_standard,
        'vert-mito'   => \&_table_vert_mito,
        'inverte-mito'=> \&_table_inverte_mito,
        'plant-plastid'=> \&_table_plant_plastid,
    );

    unless (exists $tables{$name}) {
        die "Error: unknown codon table '$name'. Supported: "
          . join(', ', sort keys %tables) . "\n";
    }
    return $tables{$name}->();
}

# NCBI Translation Table 1: Standard
sub _table_standard {
    my @base = qw(A C G T);
    my %aa = (
        # Phe
        TTT => 'F', TTC => 'F', TTA => 'L', TTG => 'L',
        # Leu
        CTT => 'L', CTC => 'L', CTA => 'L', CTG => 'L',
        # Ile / Met
        ATT => 'I', ATC => 'I', ATA => 'I', ATG => 'M',
        # Val
        GTT => 'V', GTC => 'V', GTA => 'V', GTG => 'V',
        # Ser
        TCT => 'S', TCC => 'S', TCA => 'S', TCG => 'S',
        # Pro
        CCT => 'P', CCC => 'P', CCA => 'P', CCG => 'P',
        # Thr
        ACT => 'T', ACC => 'T', ACA => 'T', ACG => 'T',
        # Ala
        GCT => 'A', GCC => 'A', GCA => 'A', GCG => 'A',
        # Tyr / Stop
        TAT => 'Y', TAC => 'Y', TAA => '*', TAG => '*',
        # His / Gln
        CAT => 'H', CAC => 'H', CAA => 'Q', CAG => 'Q',
        # Asn / Lys
        AAT => 'N', AAC => 'N', AAA => 'K', AAG => 'K',
        # Asp / Glu
        GAT => 'D', GAC => 'D', GAA => 'E', GAG => 'E',
        # Cys / Stop / Trp
        TGT => 'C', TGC => 'C', TGA => '*', TGG => 'W',
        # Arg
        CGT => 'R', CGC => 'R', CGA => 'R', CGG => 'R',
        # Ser (CTA/CTG already Leu in standard)
        AGT => 'S', AGC => 'S',
        # Arg (AGA/AGG)
        AGA => 'R', AGG => 'R',
        # Gly
        GGT => 'G', GGC => 'G', GGA => 'G', GGG => 'G',
    );
    return \%aa;
}

# NCBI Translation Table 2: Vertebrate Mitochondrial
sub _table_vert_mito {
    my $t = _table_standard();
    # Differences from standard
    $t->{AGA} = '*';  # Arg → Stop
    $t->{AGG} = '*';  # Arg → Stop
    $t->{ATA} = 'M';  # Ile → Met
    $t->{TGA} = 'W';  # Stop → Trp
    return $t;
}

# NCBI Translation Table 5: Invertebrate Mitochondrial
sub _table_inverte_mito {
    my $t = _table_standard();
    # Differences from standard
    $t->{AGA} = 'S';  # Arg → Ser
    $t->{AGG} = 'S';  # Arg → Ser
    $t->{ATA} = 'M';  # Ile → Met
    $t->{TGA} = 'W';  # Stop → Trp
    return $t;
}

# NCBI Translation Table 11: Bacterial/Archaeal/Plant Plastid
# (Identical to Standard for the coding region differences)
sub _table_plant_plastid {
    my $t = _table_standard();
    # Only differs in alternative initiation codons, which don't affect 4DTv
    return $t;
}
