=pod

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

=head1 NAME

    Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetFiles

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

