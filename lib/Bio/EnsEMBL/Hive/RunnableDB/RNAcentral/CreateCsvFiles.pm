
=pod

=head1 NAME



=head1 SYNOPSIS


=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::CreateCsvFiles;

use strict;

use Bio::RNAcentral::InputFiles;
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

    # get ncr file
    my $input_file = $self->param_required('ncr_file');

    my $opt = {};
    $opt->{'out'}      = $self->param_required('out');
    $opt->{'user'}     = $self->param_required('oracle-user');
    $opt->{'password'} = $self->param_required('oracle-password');
    $opt->{'sid'}      = $self->param_required('oracle-sid');
    $opt->{'port'}     = $self->param_required('oracle-port');
    $opt->{'host'}     = $self->param_required('oracle-host');

    # produce csv files
    my $rnac = Bio::RNAcentral::InputFiles->new($opt);
    my @files = $rnac->embl2csv($input_file);

    # my @csv_files = map { { 'csv_file' => $_ } } values @files;

        # store them for future use:
    # $self->param('csv_files', \@csv_files);

    # # import csv files
    # my $sqlldr = Bio::RNAcentral::SqlldrImport->new($opt);
    # for my $csv_file (@csv_files) {
    #     $sqlldr->load_seq($csv_file);
    # }
}

=head2 write_output

    Description :

=cut

sub write_output {
    # my $self = shift @_;

    # my $files = $self->param('csv_files');

    # $self->dataflow_output_id($files, 1);
}


1;

