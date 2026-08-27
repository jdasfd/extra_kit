#!/usr/bin/perl -w
#
# xlsx2table.pl -- Read xlsx/xls content and export as csv/tsv/md tables
#
# Unified from: xlsx2csv.pl (single file -> csv), xlsx2xls.pl (batch converting)
# and xlsx_table.pl (cell range extraction). The Win32::OLE based .xls writing
# is dropped, since this script focuses on reading content for downstream
# processing.
#
# Author: Yuqian Jiang
# Created: 2026-08-27

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Path::Tiny;
use File::Find::Rule;
use Spreadsheet::XLSX;
use Spreadsheet::ParseExcel;
use Text::CSV_XS;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

xlsx2table.pl - read xlsx/xls content and export as csv/tsv/md

=head1 SYNOPSIS

    perl xlsx2table.pl [options]
      Options:
        --help          -h          brief help message
        --file          -f  STR     one xlsx/xls file (single-file mode)
        --dir           -d  STR     dir containing excel files (batch mode)
        --suffix        -s  STR     file glob in batch mode, default is [*.xlsx]
        --sheet             STR     one sheet name, default is the first sheet
        --range             STR     cell range to extract, e.g. A1:C10
        --format        -t  STR     output format: csv|tsv|md, default is [csv]
        --output        -o  STR     output filename (single-file mode) or output
                                    dir (batch mode); default is stdout in
                                    single-file mode, or beside inputs in batch
                                    mode

    perl xlsx2table.pl -f stat.xlsx
    perl xlsx2table.pl -f stat.xlsx -t md -o stat.md
    perl xlsx2table.pl -f stat.xlsx --sheet Sheet2 --range A1:D20 -t tsv
    perl xlsx2table.pl -d ./charts -s "*.xlsx" -t tsv -o ./tables

    Notes:
        .xlsx files are parsed by Spreadsheet::XLSX, .xls files by
        Spreadsheet::ParseExcel. Cell values are raw (unformatted) values.
        In md format, the first row of the extracted table is the header.

=cut

GetOptions(
    'help|h'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \( my $file_excel ),
    'dir|d=s'    => \( my $dir ),
    'suffix|s=s' => \( my $suffix = "*.xlsx" ),
    'sheet=s'    => \( my $sheetname ),
    'range=s'    => \( my $range ),
    'format|t=s' => \( my $format = 'csv' ),
    'output|o=s' => \( my $output ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
binmode STDOUT, ':encoding(UTF-8)';

$format = lc $format;
if ( $format !~ /^(?:csv|tsv|md)$/ ) {
    die "Error: invalid format [$format], choose one of csv|tsv|md\n";
}

if ( defined $file_excel && defined $dir ) {
    die "Error: use either --file or --dir, not both\n";
}
if ( !defined $file_excel && !defined $dir ) {
    die Getopt::Long::HelpMessage(1);
}

my @files;
if ( defined $file_excel ) {
    if ( !path($file_excel)->is_file ) {
        die "Error: can't find file [$file_excel]";
    }
    @files = ( path($file_excel)->absolute->stringify );
}
else {
    if ( !path($dir)->is_dir ) {
        die "Error: can't find dir [$dir]";
    }
    @files = map { path($_)->absolute->stringify }
        File::Find::Rule->file->name($suffix)->in($dir);
    printf "\n----Total %s Files: %4s----\n\n", $suffix, scalar @files;

    if ( defined $output ) {
        path($output)->mkpath if !path($output)->is_dir;
    }
}

#----------------------------------------------------------#
# Main program
#----------------------------------------------------------#
for my $file (@files) {
    my $workbook = load_workbook($file);
    my @sheets   = $workbook->worksheets;
    if ( !@sheets ) {
        warn "No sheets in [$file], skipped\n";
        next;
    }

    my $sheet;
    if ( defined $sheetname ) {
        ($sheet) = grep { $_->get_name eq $sheetname } @sheets;
        if ( !defined $sheet ) {
            warn "Can't find sheet [$sheetname] in [$file], skipped\n";
            next;
        }
    }
    else {
        $sheet = $sheets[0];
    }

    my @rows = extract_rows( $sheet, $range );

    # decide the output target
    my ( $out_fh, $outfile );
    if ( defined $dir ) {    # batch mode, one output file per input
        my $base = path($file)->basename;
        $base =~ s/\.[^.]+$//;
        $outfile
            = defined $output
            ? path( $output, "$base.$format" )->stringify
            : path( path($file)->parent, "$base.$format" )->stringify;
    }
    elsif ( defined $output ) {    # single-file mode with explicit output
        $outfile = $output;
    }

    if ( defined $outfile ) {
        print "Write: $outfile\n";
        $out_fh = path($outfile)->filehandle( '>', ':encoding(UTF-8)' );
    }
    else {
        $out_fh = *STDOUT;
    }

    if ( !@rows ) {
        warn "Sheet [" . $sheet->get_name . "] in [$file] is empty\n";
    }
    elsif ( $format eq 'md' ) {
        write_md( $out_fh, \@rows );
    }
    else {
        write_csv( $out_fh, \@rows, $format eq 'tsv' ? "\t" : "," );
    }

    if ( defined $outfile ) {
        close $out_fh;
    }
}

exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub load_workbook {
    my $file = shift;

    if ( $file =~ /\.xlsx$/i ) {
        return Spreadsheet::XLSX->new($file);
    }
    else {
        return Spreadsheet::ParseExcel->new->parse($file);
    }
}

sub extract_rows {
    my $sheet = shift;
    my $range = shift;

    my ( $row_min, $col_min, $row_max, $col_max );
    if ( defined $range ) {
        ( $row_min, $col_min, $row_max, $col_max ) = range_to_rowcol($range);
    }
    else {
        return () unless defined $sheet->{MaxRow};
        $row_min = $sheet->{MinRow};
        $col_min = $sheet->{MinCol};
        $row_max = $sheet->{MaxRow};
        $col_max = $sheet->{MaxCol} // $sheet->{MinCol};
    }

    my @rows;
    for my $row ( $row_min .. $row_max ) {
        my @fields;
        for my $col ( $col_min .. $col_max ) {
            my $cell = $sheet->{Cells}[$row][$col];
            push @fields, defined $cell ? $cell->{Val} : q{};
        }
        push @rows, \@fields;
    }

    return @rows;
}

sub write_csv {
    my $fh   = shift;
    my $rows = shift;
    my $sep  = shift;

    my %opt = ( binary => 1, eol => "\n", sep_char => $sep );
    if ( $sep eq "\t" ) { $opt{quote_char} = undef; }   # plain tsv, no quotes
    my $csv = Text::CSV_XS->new( \%opt );

    for my $row ( @{$rows} ) {
        $csv->say( $fh, $row );
    }
}

sub write_md {
    my $fh   = shift;
    my $rows = shift;

    # pad all rows to the same column number
    my $ncol = 0;
    for my $row ( @{$rows} ) {
        $ncol = scalar @{$row} if scalar @{$row} > $ncol;
    }
    for my $row ( @{$rows} ) {
        push @{$row}, q{} while scalar @{$row} < $ncol;
    }

    my $header = shift @{$rows};
    print $fh "| ", join( " | ", map { md_escape($_) } @{$header} ), " |\n";
    print $fh "| ", join( " | ", ("---") x $ncol ), " |\n";
    for my $row ( @{$rows} ) {
        print $fh "| ", join( " | ", map { md_escape($_) } @{$row} ), " |\n";
    }
}

sub md_escape {
    my $text = shift;

    $text = q{} unless defined $text;
    $text =~ s/\r?\n/ /g;    # newlines break markdown tables
    $text =~ s/\|/\\|/g;

    return $text;
}

# "A1:C10" or "B5" -> ( row1, col1, row2, col2 ), zero-indexed
sub range_to_rowcol {
    my $range = shift;

    my @cells = split /:/, $range;
    if ( @cells == 1 ) {
        return ( cell_to_rowcol( $cells[0] ), cell_to_rowcol( $cells[0] ) );
    }
    elsif ( @cells == 2 ) {
        return ( cell_to_rowcol( $cells[0] ), cell_to_rowcol( $cells[1] ) );
    }
    else {
        die "Error: invalid range [$range], expected like A1:C10\n";
    }
}

# "C10" -> ( row, col ), zero-indexed
sub cell_to_rowcol {
    my $cell = shift;

    if ( !defined $cell || $cell !~ /^([A-Za-z]{1,3})(\d+)$/ ) {
        die "Error: invalid cell [$cell], expected like A1\n";
    }

    my $col_str = uc $1;
    my $row     = $2;

    my $col = 0;
    for my $char ( split //, $col_str ) {    # base26 string to number
        $col = $col * 26 + ( ord($char) - ord('A') + 1 );
    }

    return ( $row - 1, $col - 1 );
}

__END__
