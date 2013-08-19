
=pod

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::ImportCsv;

use strict;

use Bio::RNAcentral::SqlldrImport;

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

    my $input_file = $self->param_required('csv_file');

    my $opt = {};
    $opt->{'output_folder'} = $self->param_required('output_folder');
    $opt->{'user'}     = $self->param_required('oracle-user');
    $opt->{'password'} = $self->param_required('oracle-password');
    $opt->{'sid'}      = $self->param_required('oracle-sid');
    $opt->{'port'}     = $self->param_required('oracle-port');
    $opt->{'host'}     = $self->param_required('oracle-host');

    my $sqlldr = Bio::RNAcentral::SqlldrImport->new($opt);
    $sqlldr->load_seq($input_file);
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
}

1;

