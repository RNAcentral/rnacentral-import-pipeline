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

=head1 DESCRIPTION

    Recursively get a list of files with a specified extension.

=cut

package Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetNcProduct;

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

    my $opt = {};
    $opt->{'output_folder'} = $self->param_required('output_folder');
    $opt->{'release_type'}  = $self->param_required('release_type');

    my $rnac = Bio::RNAcentral::InputFiles->new($opt);
    my @files = $rnac->list_folder_recursive_ftp();

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

