=pod

=head1 NAME

    Bio::RNAcentral::InputFiles

=head1 SYNOPSIS

    my $opt{'out'} = 'path/to/output/temp/directory';
    my $a = Bio::RNAcentral::InputFiles->new($opt);
    # prepare csv files for sqlldr
    $a->process_folder('path/to/files1');
    $a->process_folder('path/to/files2');

=head1 DESCRIPTION

    This package is responsible for getting a list of input files,
    reorganizing them as necessary, and parsing them into csv files
    expected by SQL*Loader.

=head1 CONTACT


=cut

package Bio::RNAcentral::InputFiles;

use strict;
use warnings;

# override default modules, use embl.pm and CRC64.pm modules from the lib directory
use Cwd            ();
use File::Basename ();
BEGIN {
    my $root_dir = File::Basename::dirname( File::Basename::dirname( Cwd::realpath($0) ) );
    unshift @INC, "$root_dir/lib";
}

use File::Spec;
use File::Find qw(finddepth);

use Bio::SeqIO;   # BioPerl is used for reading embl files
use SWISS::CRC64; # cyclic redundancy check
use Digest::MD5 qw(md5_hex);

use base ('Bio::RNAcentral::Base');


sub new {
    my ($class, $opt) = @_;

    # run parent constructor
    my $self = $class->SUPER::new();

    # location of temporary output files
    $self->{'temp_dir'} = $opt->{'out'};

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


=head2 list_folder_recursive

    Find all files with the specified extension in a folder recursively.

=cut

sub list_folder_recursive {
    (my $self, my $dir, my $extension) = @_;

    my @files;

    finddepth(sub {
      return if(-d $_ || $_ !~ /$extension$/ || $_ =~ m/^\./ );
      push @files, $File::Find::name;
    }, $dir);

    $self->{'logger'}->info("Found " . scalar(@files) . " $extension files");

    return @files;
}


=head2 get_files

    Read all files with specified extension in a location.
    Produce csv files for sqlldr import.

=cut

sub get_files {
    (my $self, my $location) = @_;

    my @files = ();

    # get a list of all files
    my @original_files = $self->list_folder_recursive($location, $self->{'opt'}{'file_extension'});

    # split large files into small chunks and analyze them instead
    foreach my $file (@original_files) {
        @files = (@files, $self->file2chunks($file));
    }

    # @files = @original_files;

    return @files;
}


=head2 list_folder

    Get all files with the specified extension from the input directory.

=cut

sub list_folder {

    (my $self, my $dir, my $extension) = @_;

    opendir(DIR, $dir) or die $!;

    my @files = ();
    while (my $file = readdir(DIR)) {
        next if ($file =~ m/^\./);
        next if ($file !~ m/\.$extension$/);

        push @files, File::Spec->catfile($dir, $file);
    }

    closedir(DIR);

    $self->{'logger'}->info("Found " . scalar(@files) . " .$extension files");

    return @files;
}


=head2 list_folder_order_by_size

    Get all files with the specified extension from the input directory.
    In addition, return a hash with filesizes (in bytes).
    Used in conjunction with group_files.

=cut

sub list_folder_order_by_size {

    my ($self, $dir, $extension) = @_;

    opendir(DIR, $dir) or die $!;

    my %size = ();

    while (my $file = readdir(DIR)) {
        next if ($file =~ m/^\./);
        next if ($file !~ m/\.$extension$/);

        $size{$file} = -s $dir . '/' . $file;
    }
    closedir(DIR);

    my @ordered_files = ();
    foreach my $key (sort { $size{$a} <=> $size{$b} } keys %size) {
        push @ordered_files, $key;
    }

    $self->{'logger'}->info("Found " . scalar(keys %size) . " .$extension files");

    return (\%size, \@ordered_files);
}


=head2 group_files

    Return an array of arrays with each nested array representing
    a group of files with cumulative size not exceeding
    the specified cutoff.
    Called recursively. Used in conjunction with sub list_folder_order_by_size.

=cut

sub group_files {

    my $self          = shift;
    my %size          = %{$_[0]};
    my @ordered_files = @{$_[1]};

    # no files left
    if ( !scalar(@ordered_files) ) {
        return [];
    }

    my $MAX = $self->{'opt'}{'file_size_cutoff'};

    my $i = 0;
    my $file1 = shift @ordered_files;
    my @group = ( $file1 );
    my $total = $size{$file1};

    foreach my $file2 (@ordered_files) {
        if ( $size{$file2} + $total <= $MAX ) {
            push @group, $file2;
            $total += $size{$file2};
        } else {
            last; # don't check the remaining files because the array is sorted by size
        }
        $i++;
    }

    # get the remaining non-concatenated files
    my @next_level = @ordered_files[$i .. $#ordered_files];

    # recursive call
    my @x = @{$self->group_files(\%size, \@next_level)};

    # make sure the results are on the same array depth
    my @temp = (\@group); # add the group from this iteration
    foreach my $xx (@x) {
        if ( $xx ) { # skip undef from the last level
            push @temp, $xx;
        }
    }

    # contains this group and all groups from following iterations
    return \@temp;
}


=head2 file2chunks

    Split large files into chunks for faster processing.
    Return an array of filenames or an array with the original file if unmodified.

=cut

sub file2chunks {

    my $self = shift;
    my $filename = shift;

    my $size = -s $filename;
    if ( $size < $self->{'opt'}{'file_size_cutoff'} ) {
        return ($filename);
    }

    $self->{'logger'}->info("File $filename is $size bytes, will split it into chunks");

    local $/ = "//\n";

    my $i      = 1;
    my $text   = '';
    my @chunks = ();
    my $chunkfile  = '';

    my $file = File::Basename::fileparse($filename, qr/\.[^.]*/);

    use bytes; # to calculate length in bytes

    open (INFILE, $filename);
    while (<INFILE>) {
        $text .= $_;
        if ( length($text) > $self->{'opt'}{'file_size_cutoff'} ) {
            # create chunk files in temp directory
            $chunkfile = File::Spec->catfile($self->{'temp_dir'}, $file . '_chunk' . $i . '.' . $self->{'opt'}{'file_extension'});
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
        $chunkfile = File::Spec->catfile($self->{'temp_dir'}, $file . '_chunk' . $i . '.' . $self->{'opt'}{'file_extension'});
        push @chunks, $chunkfile;
        _print_to_file($chunkfile, $text);
        $self->{'logger'}->info("Created file $chunkfile");
    }

    return @chunks;
}


sub _print_to_file {
    open (OUTFILE, "> $_[0]") or die("Couldn't open file $_[0]");
    print OUTFILE $_[1];
    close OUTFILE;
}


=head2 embl2csv

    Read an embl file and create two csv files,
    one with short, one with long sequences.
    In addition, create csv files with database cross-references
    and literature references.

=cut

sub embl2csv {
    (my $self, my $file) = @_;

    # get file name without extension
    my $job_id = File::Basename::fileparse($file, qr/\.[^.]*/);

    my $fname_long  = $self->get_output_filename($self->{'temp_dir'}, $job_id, 'long',  'csv');
    my $fname_short = $self->get_output_filename($self->{'temp_dir'}, $job_id, 'short', 'csv');
    my $fname_xref  = $self->get_output_filename($self->{'temp_dir'}, $job_id, 'xref',  'csv');
    my $fname_refs  = $self->get_output_filename($self->{'temp_dir'}, $job_id, 'refs',  'csv');

    # open output files
    my $fh_long  = IO::File->new("> $fname_long")  or $self->{'logger'}->logdie("Couldn't open file $fname_long");
    my $fh_short = IO::File->new("> $fname_short") or $self->{'logger'}->logdie("Couldn't open file $fname_short");
    my $fh_xrefs = IO::File->new("> $fname_xref")  or $self->{'logger'}->logdie("Couldn't open file $fname_xref");
    my $fh_refs  = IO::File->new("> $fname_refs")  or $self->{'logger'}->logdie("Couldn't open file $fname_refs");

    # counters
    my ($i, $records, $seqs_long, $seqs_short, $dblinks_num, $refs_num);
    $i = $records = $seqs_long = $seqs_short = $dblinks_num = $refs_num = 0;

    # open file with BioPerl
    $self->{'logger'}->info("Reading $file");
    my $stream = Bio::SeqIO->new(-file => $file, -format => 'EMBL');

    my ($data, $seq, $dblinks, $refs);
    $dblinks = $refs = '';

    # loop over records
    while ( eval { $seq = $stream->next_seq() } ) {
        $i++;

        # get basic data for UniParc-style functionality
        $data = $self->_get_basic_data($seq);

        unless ( $data->{'isValid'} ) {
            $self->{'logger'}->logwarn("Skipping record $i: $data->{'text'}");
            next;
        }

        # print the data in different files depending on sequence length
        if ( $data->{'isLong'} ) {
            print $fh_long $data->{'text'};
            $seqs_long += $data->{'seqs_long'};
        } else {
            print $fh_short $data->{'text'};
            $seqs_short += $data->{'seqs_short'};
        }

        # get database links
        $dblinks = $self->_get_dblinks($seq);
        if ( $dblinks ) {
            print $fh_xrefs $dblinks->{'text'};
            $dblinks_num += $dblinks->{'num'};
        }

        # get literature references
        $refs = $self->_get_references($seq);
        if ( $refs ) {
            print $fh_refs $refs->{'text'};
            $refs_num += $refs->{'num'};
        }
    }

    # report any problems, delete already created files to prevent any import
    if ( $@ ) {
        $self->{'logger'}->logwarn('BioPerl error');
        $self->{'logger'}->logwarn($@);
        $self->{'logger'}->logwarn("Not safe to continue, skipping $file");
        unlink $fname_short, $fname_long, $fname_refs, $fname_xref;
        return ();
    }

    $records += $i;

    $self->{'logger'}->info("$records records, $seqs_short short, $seqs_long long");
    $self->{'logger'}->info("Created files $fname_long and $fname_short");
    $self->{'logger'}->info("Found $dblinks_num database links");
    $self->{'logger'}->info("Found $refs_num literature references");

    $fh_short->close;
    $fh_long->close;
    $fh_refs->close;
    $fh_xrefs->close;

    $self->_check_temp_files($fname_xref);
    $self->_check_temp_files($fname_refs);

    # return an array with csv files
    return ($self->_check_temp_files($fname_short), $self->_check_temp_files($fname_long));
}


# delete temp files if empty, return only files with content
sub _check_temp_files{
    (my $self, my $file) = @_;

    if ( -s $file > 0 ) {
        return $file;
    } else {
        $self->{'logger'}->info("File $file is empty and will be deleted");
        unlink $file;
        return ();
    }
}


=head2 _get_references

    Get csv text with literature references given a BioPerl seq object.

=cut

sub _get_references {

    (my $self, my $seq)  = @_;

    my ($ref_authors, $ref_location, $ref_title,
        $ref_pubmed,  $ref_doi,      $ref_publisher,
        $ref_editors, $ref_consortium);
    my $text = '';
    my $num = 0;

    my $anno_collection = $seq->annotation;
    for my $key ( $anno_collection->get_all_annotation_keys ) {
        my @annotations = $anno_collection->get_Annotations($key);
        for my $value ( @annotations ) {
            # all reference lines: RN, RC, RP, RX, RG, RA, RT, RL
            if ( $value->tagname eq 'reference' ) {
                $num++;
                $ref_authors   = _nvl($value->authors());   # RA line
                $ref_location  = _nvl($value->location());  # RL line
                $ref_title     = _nvl($value->title());     # RT line
                $ref_pubmed    = _nvl($value->pubmed());    # parsed RX line
                $ref_doi       = _nvl($value->doi());       # parsed RX line
                $ref_publisher = _nvl($value->publisher()); # parsed RL line
                $ref_editors   = _nvl($value->editors());   # parsed RL line
                $text .= '"' . join('","', ($seq->display_id,
                                            $ref_authors,
                                            $ref_location,
                                            $ref_title,
                                            $ref_pubmed,
                                            $ref_doi,
                                            $ref_publisher,
                                            $ref_editors) ) . "\"\n";
            }
        }
    }
    return { text => $text, num => $num };
}


=head2 _get_dblinks

    Get csv text with database crossreferences given a BioPerl seq object.

=cut

sub _get_dblinks {

    (my $self, my $seq)  = @_;

    my ($database, $primary_id, $optional_id);
    my $text = '';
    my $num = 0;

    # reading database links (DR lines) and references
    my $anno_collection = $seq->annotation;
    for my $key ( $anno_collection->get_all_annotation_keys ) {
        my @annotations = $anno_collection->get_Annotations($key);
        for my $value ( @annotations ) {
            if ( $value->tagname eq "dblink" ) {
                $num++;
                $database    = _nvl($value->database());
                $primary_id  = _nvl($value->primary_id());
                $optional_id = _nvl($value->optional_id());
                $text .= '"' . join('","', ($seq->display_id,
                                            $database,
                                            $primary_id,
                                            $optional_id) ) . "\"\n";
            }
        }
    }

    return { text => $text, num => $num };
}


=head2 _get_references

    Get minimum data required for UniParc-style functionality given a BioPerl seq object.

=cut

sub _get_basic_data {

    (my $self, my $seq)  = @_;

    my ($ac, $length, $version, $isLong, $md5, $crc64, $taxid,
        $data, $sequence, $seqs_short, $seqs_long, $isValid);

    $ac = $md5 = $crc64 = $taxid = $data = $sequence = '';
    $length = $version = $seqs_long = $seqs_short = $isValid = 0;

    # get new data
    $ac       = $seq->display_id;
    $length   = $seq->length;
    $version  = $seq->seq_version;
    $sequence = $seq->seq;

    if ($length > $self->{'opt'}{'maxseqshort'}) {
        $isLong = 1;
        $seqs_long++;
    } else {
        $isLong = 0;
        $seqs_short++;
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

    # quality control
    if ( $ac eq '' or $sequence eq '' or $length == 0 ) {
        $isValid = 0;
    } else {
        $isValid = 1;
    }

    return {
        'text'   => join(',', ($crc64, $length, $sequence, $ac, $version, $taxid, $md5)) . "\n",
        'isLong' => $isLong,
        'seqs_long'  => $seqs_long,
        'seqs_short' => $seqs_short,
        'isValid'    => $isValid,
    };
}

# get non-empty value, named after the NVL function in Oracle
sub _nvl {
    my $value = shift;
    return ( defined( $value ) ? $value : '' );
}


1;