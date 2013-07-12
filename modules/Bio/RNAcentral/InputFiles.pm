#!/usr/bin/env perl

package Bio::RNAcentral::InputFiles;

=pod

=head1 NAME

test

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT


=cut

use strict;
use warnings;

# override default modules, use embl.pm and CRC64.pm modules from the lib directory
use Cwd            ();
use File::Basename ();
BEGIN {
    my $root_dir = File::Basename::dirname( File::Basename::dirname( Cwd::realpath($0) ) );
    unshift @INC, "$root_dir/lib";
}

use Bio::SeqIO;   # BioPerl is used for reading embl files
use SWISS::CRC64; # cyclic redundancy check
use Digest::MD5 qw(md5_hex);

use base ('Bio::RNAcentral::Base');

# Global parameters
# our $MAXSEQLONG    = 1000000;  # maximum length for long sequences stored as clobs
our $MAXSEQSHORT   = 4000;     # maximum length for short sequences stored as chars
our $EXTENSION     = 'ncr';    # look for .ncr files
our $SIZECUTOFF    = 10**7;    # 10 Mb, file size cutoff


sub new {
    my ($class, @args) = @_;

    # run parent constructor
    my $self = $class->SUPER::new();

    return $self;
}


=head2 process_folder

    Read all files with specified extension in a location.
    Produce csv files for sqlldr import.

=cut

sub process_folder {
    (my $self, my $location) = @_;
    
    # concatenate small files, split large ones
    my @files = $self->get_files($location);

    # create csv files based on embl files
    for my $file (@files) {
        # TODO get a list of newly created output files
        $self->embl2csv($file);
    }
}


# single job to get files, then spawn jobs on individual files
sub get_files {
    (my $self, my $location) = @_;    

    my @files = ();

    # get a list of all files
    my @original_files = $self->list_files($location, $EXTENSION);

    # split large files into small chunks and analyze them instead
    foreach my $file (@original_files) {
        @files = (@files, $self->file2chunks($file));
    }

    #reorganize_files;

    return @files;
}


sub reorganize_files {
    
    #return @files; # similarly-sized files
}



=head2

    Read an embl file and create two csv files,
    one with short, one with long sequences.

=cut

sub embl2csv {
    (my $self, my $file, my $job_id) = @_;

    # get file name without extension
    $job_id = File::Basename::fileparse($file, qr/\.[^.]*/);

    my $fname_long  = $self->get_output_filename($job_id, 'long',  'csv');
    my $fname_short = $self->get_output_filename($job_id, 'short', 'csv');

    # open output files
    my $fh_long  = IO::File->new("> $fname_long")  or $self->logdie("Couldn't open file $fname_long");
    my $fh_short = IO::File->new("> $fname_short") or $self->logdie("Couldn't open file $fname_short");
        # get other filehandles, if necessary

    # counters
    my $i = 0;    
    my $records = 0;

    # open file with BioPerl        
    $self->{'logger'}->info("Reading $file");
    my $stream = Bio::SeqIO->new(-file => $file, -format => 'EMBL');

    my %data = ();

    # loop over records
    while ( (my $seq = $stream->next_seq()) ) {
        $i++;

        # get basic data for UniParc-style functionality
        %data = $self->_get_basic_data($seq);

        # print the data in different files depending on sequence length
        if ( $data{'isLong'} == 1 ) {
            print $fh_long $data{'text'};
        } elsif ( $data{'isLong'} == 0 ) {
            print $fh_short $data{'text'};
        } else {
            print 'Error';
        }

        # get other data from the same record 
        # and print it to other files if necessary
    }
    $records += $i;
    
    $self->{'logger'}->info("$records records, $data{'seqs_short'} short, $data{'seqs_long'} long");
    $self->{'logger'}->info("Created file $fname_long and $fname_short");

    $fh_short->close;
    $fh_long->close;

    # TODO delete long or short file if empty
}



# minimum data required for UniParc-style functionality
sub _get_basic_data {

    (my $self, my $seq)  = @_;

    my $ac       = '';
    my $length   = 0;
    my $version  = 0;
    my $md5      = '';
    my $crc64    = '';
    my $taxid    = '';
    my $isLong   = -1; # directs sequences to different files
    my $data     = '';
    my $sequence = '';

    my $seqs_long  = 0;
    my $seqs_short = 0;    

    # get new data
    $ac       = $seq->display_id;
    $length   = $seq->length;
    $version  = $seq->seq_version;
    $sequence = $seq->seq;

    if ($length > $MAXSEQSHORT) {
        $isLong = 1;
        $seqs_long++;
    } else {
        $seqs_short++;
        $isLong = 0;
    }

    $md5   = md5_hex($seq->seq);
    $crc64 = SWISS::CRC64::crc64($seq->seq);

    for my $feat_object ($seq->get_SeqFeatures) {
        if ( $feat_object->primary_tag eq 'source' ) {
            for my $tag ($feat_object->get_all_tags) {
                for my $value ($feat_object->get_tag_values($tag)) {
                    if ($tag eq 'db_xref') {
                        (my $text, $taxid) = split(':', $value);
                        last;
                    }
                }
            }
            last; # found taxid, stop looking
        }
    }

    return (
        'text'   => join(',', ($crc64, $length, $sequence, $ac, $version, $taxid, $md5)) . "\n",
        'isLong' => $isLong,
        'seqs_long'  => $seqs_long,
        'seqs_short' => $seqs_short,
    );
}


=head2 list_files

    Get all files with the specified extension from the input directory

=cut

sub list_files {

    (my $self, my $dir, my $extension) = @_;

    opendir(DIR, $dir) or die $!;

    my @files = ();
    while (my $file = readdir(DIR)) {
        next if ($file =~ m/^\./);
        next if ($file !~ m/\.$extension$/);

        push @files, $dir . '/' . $file;
    }

    closedir(DIR);

    $self->{'logger'}->info("Found " . scalar(@files) . " .$extension files");

    return @files;
}


=head2 file2chunks

    Some input files are very large, so it's faster to split them into chunks for processing.
    Split large files into chunks and return an array of filenames.
    Return an array with the original file if it's small.

=cut

sub file2chunks {

    my $self = shift;
    my $filename = shift;

    my $size = -s $filename;
    if ( $size < $SIZECUTOFF ) {
        return ($filename);
    }

    $self->{'logger'}->info("File $filename is $size bytes, will split it into chunks");

    local $/ = "//\n";

    my $i      = 1;
    my $text   = '';
    my @chunks = ();
    my $chunkfile  = '';
    (my $path, my $extension) = split('\.', $filename);

    use bytes; # to calculate length in bytes

    open (INFILE, $filename);
    while (<INFILE>) {
        $text .= $_;
        if ( length($text) > $SIZECUTOFF ) {
            $chunkfile = $path . '_chunk' . $i . '.' . $extension;
            push @chunks, $chunkfile;
            _print_to_file($chunkfile, $text);
            $i++;
            $text = '';
            $self->{'logger'}->info("Created file $chunkfile");
        }
    }
    close INFILE;

    # left over sequences
    if ( length($text) > 0 ) {
        $chunkfile = $path . '_chunk' . $i . '.' . $extension;
        push @chunks, $chunkfile;
        _print_to_file($chunkfile, $text);
        $self->{'logger'}->info("Created file $chunkfile");
    }

    return @chunks;
}


sub _print_to_file {
    open (OUTFILE, "> $_[0]") or die("Couldn't open file");
    print OUTFILE $_[1];
    close OUTFILE;
}

1;