
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

