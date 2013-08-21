=pod

=head1 NAME

    Bio::RNAcentral::InputFiles

=head1 SYNOPSIS


=head1 DESCRIPTION

    This package is responsible for getting a list of input files,
    reorganizing them as necessary.

=head1 CONTACT


=cut

package Bio::RNAcentral::InputFiles;

use strict;
use warnings;

use File::Basename ();
use File::Spec;
use File::Find qw(finddepth);

use base ('Bio::RNAcentral::Base');


sub new {
    my ($class, $opt) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt);

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
            $chunkfile = File::Spec->catfile($self->{'output_folder'}, $file . '_chunk' . $i . '.' . $self->{'opt'}{'file_extension'});
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
        $chunkfile = File::Spec->catfile($self->{'output_folder'}, $file . '_chunk' . $i . '.' . $self->{'opt'}{'file_extension'});
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


1;