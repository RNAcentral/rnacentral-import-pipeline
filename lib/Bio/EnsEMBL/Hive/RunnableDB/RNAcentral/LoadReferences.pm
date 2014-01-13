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

    Load literature references extracted from the Non-coding Product.

=cut

package Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::LoadReferences;

use strict;

use Bio::RNAcentral::SqlldrImportReferences;

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
}


=head2 run

    Description :


=cut

sub run {
    my $self = shift @_;

    my $opt = {};
    $opt->{'user'}          = $self->param_required('oracle-user');
    $opt->{'password'}      = $self->param_required('oracle-password');
    $opt->{'sid'}           = $self->param_required('oracle-sid');
    $opt->{'port'}          = $self->param_required('oracle-port');
    $opt->{'host'}          = $self->param_required('oracle-host');
    $opt->{'output_folder'} = $self->param_required('output_folder');

    my $rnac = Bio::RNAcentral::SqlldrImportReferences->new($opt, 'refs');
    $rnac->update();
}


=head2 write_output

    Description :


=cut

sub write_output {
}

1;

