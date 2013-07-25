
=pod

=head1 NAME



=head1 SYNOPSIS


=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::SqlldrUpload;

use strict;

use Bio::RNAcentral::SqlldrImport;

use base ('Bio::EnsEMBL::Hive::Process');


=head2 param_defaults

    Description : Implements param_defaults() interface method of Bio::EnsEMBL::Hive::Process that defines module defaults for parameters.

=cut

sub param_defaults {

    return {
    };
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

    # get csv file
    my $input_file = $self->param_required('csv_file');

    print $input_file , "\n\n";

    my $opt = {};
    $opt->{'out'}      = $self->param_required('out');
    $opt->{'user'}     = $self->param_required('oracle-user');
    $opt->{'password'} = $self->param_required('oracle-password');
    $opt->{'sid'}      = $self->param_required('oracle-sid');
    $opt->{'port'}     = $self->param_required('oracle-port');
    $opt->{'host'}     = $self->param_required('oracle-host');


    # produce one or two csv files
    my $sqlldr = Bio::RNAcentral::SqlldrImport->new($opt);
    $sqlldr->load_seq($input_file);

}

=head2 write_output

    Description :


=cut

sub write_output {


}


1;

