
=pod

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 CONTACT

=cut


package Bio::EnsEMBL::Hive::PipeConfig::RNAcentralUpdate_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description :


=cut

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'rnacentral_update',             # name used by the beekeeper to prefix job names on the farm
        'hive_force_init' => 1,                             # always recreate the hive database

        # get the command line option
        'location'        => $self->o('in'),
    };
}


=head2 pipeline_create_commands

    Description :


=cut

sub pipeline_create_commands {
    my ($self) = @_;
    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation

    ];
}


=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.

=cut

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class

        # store command line parameters
        'out'             => $self->o('out'),
        'oracle-user'     => $self->o('oracle-user'),
        'oracle-password' => $self->o('oracle-password'),
        'oracle-sid'      => $self->o('oracle-sid'),
        'oracle-port'     => $self->o('oracle-port'),
        'oracle-host'     => $self->o('oracle-host'),
    };
}


=head2 pipeline_analyses

    Description :

=cut

sub pipeline_analyses {
    my ($self) = @_;

    return [
        {   -logic_name => 'get_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetFiles',
            -meadow_type=> 'LOCAL',
            -analysis_capacity  =>  1,
            -input_ids => [
                { 'location' => $self->o('location') },
            ],
            -flow_into => {
                1 => { 'print_files' => { 'ncr_file' => '#ncr_file#' } },
            },
        },

        {   -logic_name    => 'print_files',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::CreateCsvFiles',
            -analysis_capacity  =>  4,
            -flow_into => {
                1 => { 'sqlldr' => { 'csv_file' => '#csv_file#' } },
            },
        },

        {   -logic_name    => 'sqlldr',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::SqlldrUpload',
            -analysis_capacity  =>  4,
        },

    ];
}

1;

