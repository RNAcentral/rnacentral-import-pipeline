=head1 LICENSE

Copyright [2009-2014] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

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
use File::Listing qw(parse_dir);
use Net::FTP;

use base ('Bio::RNAcentral::Base');


sub new {
    my ($class, $opt) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt);

    return $self;
}


=head2

    Download the entire non-coding product and uncompress all files.
    Split large files into chunks.
    Return a list of files for further processing.

=cut

sub list_folder_recursive_ftp {
    my $self = shift;

    my @output_files = ();

    $self->{'logger'}->info("Connecting to the non-coding product FTP site");

    # create a new instance of the FTP connection
    my $ftp = Net::FTP->new($self->{'opt'}{'ebi_ftp_site'}, Debug=>0)
        or $self->{'logger'}->logdie("Cannot connect $!");

    # login to the server
    $ftp->login($self->{'opt'}{'ebi_ftp_user'}, $self->{'opt'}{'ebi_ftp_password'})
        or $self->{'logger'}->logdie("Login failed $!");

    # change remote directory
    my $remote_dir;
    if ( $self->{'release_type'} eq 'F' ) {
        $self->{'logger'}->info('Getting full release');
        $remote_dir = $self->{'opt'}{'ebi_ftp_non_coding_product_release'};
    } elsif ( $self->{'release_type'} eq 'I' ) {
        $self->{'logger'}->info('Getting incremental release');
        $remote_dir = $self->{'opt'}{'ebi_ftp_non_coding_product_update'};
    } else {
        $self->{'logger'}->logdie("Incorrect release type");
    }

    $ftp->cwd($remote_dir)
        or $self->{'logger'}->logdie("Could not change remote working directory");

    # set binary transfer mode, necessary for gz files
    $ftp->binary();

    # get recursive file listing
    $self->{'logger'}->info("Getting recursive file list");
    my @ls = $ftp->ls('-lR');

    my $ftpdir = $self->get_ftp_downloads_path();
    if ( !-e $ftpdir ) {
        mkdir $ftpdir;
    }

    my ($filename, $savename, $name, $type, $size, $mtime, $mode, $uncompressed_file, $chunks);

    # parse and loop through the directory listing
    foreach my $file (parse_dir(\@ls))
    {
        ($name, $type, $size, $mtime, $mode) = @$file;

        next if $name !~ /\.$self->{'opt'}{'file_extension'}(\.gz)*$/;

        $filename = File::Basename::fileparse($name); # get just the filename
        $savename = File::Spec->catfile($ftpdir, $filename);

        $uncompressed_file = $savename;
        $uncompressed_file =~ s/(\.gz)*$//;

        $self->{'logger'}->info("Downloading $savename");

        $ftp->get("$remote_dir/$name", $savename) or $self->{'logger'}->logdie($!);

        $self->{'logger'}->info("Unzipping");

        system("gunzip -fq $savename"); # will create $uncompressed_file

        $chunks = $self->file2chunks($uncompressed_file);

        push @output_files, @$chunks;
    }

    $self->{'logger'}->info("Ftp sync complete, processed " . scalar @output_files . " files");

    # close the FTP connection
    $ftp->quit();

    return @output_files;
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
        @files = (@files, @{ $self->file2chunks($file) });
    }

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
        return [ $filename ];
    }

    $self->{'logger'}->info("File $filename is $size bytes, will split it into chunks");

    local $/ = "//\n";

    my $i      = 1;
    my $text   = '';
    my @chunks = ();
    my $chunkfile  = '';

    # get location of the chunk files
    my $chunks_path = $self->get_chunks_path();
    if ( !-e $chunks_path ) {
        mkdir $chunks_path;
    }

    my $file = File::Basename::fileparse($filename, qr/\.[^.]*/);

    use bytes; # to calculate length in bytes

    open (INFILE, $filename);
    while (<INFILE>) {
        $text .= $_;
        if ( length($text) > $self->{'opt'}{'file_size_cutoff'} ) {
            # create chunk files in temp directory
            $chunkfile = File::Spec->catfile($chunks_path, $file . '_chunk' . $i . '.' . $self->{'opt'}{'file_extension'});
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
        $chunkfile = File::Spec->catfile($chunks_path, $file . '_chunk' . $i . '.' . $self->{'opt'}{'file_extension'});
        push @chunks, $chunkfile;
        _print_to_file($chunkfile, $text);
        $self->{'logger'}->info("Created file $chunkfile");
    }

    return \@chunks;
}


sub _print_to_file {
    open (OUTFILE, "> $_[0]") or die("Couldn't open file $_[0]");
    print OUTFILE $_[1];
    close OUTFILE;
}


1;
