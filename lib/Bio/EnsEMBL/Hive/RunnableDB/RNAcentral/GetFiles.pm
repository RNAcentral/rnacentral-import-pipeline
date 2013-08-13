
=pod

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetFiles;

use strict;

use Bio::RNAcentral::InputFiles;

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
    $opt->{'out'} = $self->param_required('out');

    my $rnac = Bio::RNAcentral::InputFiles->new($opt);
    my @files = $rnac->list_folder_recursive($location, $extension);

    my @files = map { { 'ncr_file' => $_ } } values @files;

        # store them for future use:
    $self->param('ncr_files', \@files);
}


=head2 run

    Description : TODO: add file concatenation and splitting here.


=cut

sub run {
}


=head2 write_output

    Description :


=cut

sub write_output {  # nothing to write out, but some dataflow to perform:
    my $self = shift @_;

    my $files = $self->param('ncr_files');

    $self->dataflow_output_id($files, 1);
}

1;

