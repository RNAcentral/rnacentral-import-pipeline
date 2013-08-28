
=pod

=head1 DESCRIPTION

    Prepare the Oracle staging table for loading the data.

=cut

package Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::TruncateStagingTable;

use strict;


use Bio::RNAcentral::OracleUpdate;

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

    # truncate the sequences staging table
    my $oracle = Bio::RNAcentral::OracleUpdate->new($opt);
    $oracle->db_oracle_connect();
    $oracle->truncate_staging_table();
    $oracle->db_oracle_disconnect();
}


=head2 write_output

    Description :


=cut

sub write_output {
}

1;

