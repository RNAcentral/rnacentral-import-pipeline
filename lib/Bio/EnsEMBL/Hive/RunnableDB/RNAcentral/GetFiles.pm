
=pod

=head1 DESCRIPTION

    Recursively get a list of files with a specified extension.

=cut

package Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetFiles;

use strict;

use Bio::RNAcentral::InputFiles;
use File::Spec;

use base ('Bio::EnsEMBL::Hive::Process');


=head2 param_defaults

    Description :

=cut

sub param_defaults {
}


=head2 fetch_input

    Description :

=cut

sub fetch_input {
    my $self = shift @_;

    my $location = $self->param_required('location');
    my $extension = $self->param_required('extension');

    my $opt = {};
    $opt->{'output_folder'} = $self->param_required('output_folder');

    my $rnac = Bio::RNAcentral::InputFiles->new($opt);
    my @files = $rnac->list_folder_recursive($location, $extension);

    # if no files are found, create an empty dummy file
    # to enable correct hive dataflow.
    if ( !@files ) {
        my $filename = File::Spec->catfile($location, 'empty_dummy_file.' . $extension);
        open(DUMMY, ">", $filename ) || die "Can't open $filename: $!";
        close DUMMY;
        push @files, $filename;
    }

    my @files = map { { 'ncr_file' => $_ } } values @files;

        # store the files for future use:
    $self->param('ncr_files', \@files);
}


=head2 run

    Description :


=cut

sub run {
}


=head2 write_output

    Description :


=cut

sub write_output {
    my $self = shift;

    my $files = $self->param('ncr_files');

    $self->dataflow_output_id($files, 1);
}

1;

