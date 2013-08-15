
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

        'pipeline_name' => 'rnacentral_staging',             # name used by the beekeeper to prefix job names on the farm
        'hive_force_init' => 1,                             # always recreate the hive database

        # get the command line option
        'location'        => $self->o('in'),

        'pipeline_db'   => {
            -host   => $self->o('host'),
            -port   => $self->o('port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -dbname => $self->o('pipeline_name'),  # example of a linked definition (resolved via saturation)
        },

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
        'in'              => $self->o('in'),
        'oracle-user'     => $self->o('oracle-user'),
        'oracle-password' => $self->o('oracle-password'),
        'oracle-sid'      => $self->o('oracle-sid'),
        'oracle-port'     => $self->o('oracle-port'),
        'oracle-host'     => $self->o('oracle-host'),
    };
}


=head2 pipeline_analyses

    Description : Main logic of the pipeline.

=cut

sub pipeline_analyses {
    my ($self) = @_;

    return [
        {   -logic_name => 'get_ncr_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetFiles',
            -analysis_capacity  =>  1,
            -input_ids  => [
                { 'location'  => $self->o('in'),
                  'extension' => 'ncr' }
            ],
            -flow_into => {
                1 => { 'check_chunks' => { 'ncr_file' => '#ncr_file#' } },
            },
        },

        {   -logic_name => 'check_chunks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::CheckChunks',
            -flow_into => {
                1 => { 'create_csv_files' => { 'ncr_file' => '#ncr_file#' } },
            },
        },

        {   -logic_name    => 'create_csv_files',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::CreateCsvFiles',
        },

        {   -logic_name => 'get_long_csv_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetFiles',
            -analysis_capacity  => 1,
            -input_ids  => [
                { 'location'  => $self->o('out'),
                  'extension' => '_long.csv' }
            ],
            -flow_into => {
                1 => { 'import_long_csv' => { 'csv_file' => '#ncr_file#' } },
            },
            -wait_for => [ 'create_csv_files' ]
        },

        {   -logic_name    => 'import_long_csv',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::ImportCsv',
            -analysis_capacity  =>  1,
        },

        {   -logic_name => 'get_short_csv_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetFiles',
            -analysis_capacity  => 1,
            -input_ids  => [
                { 'location'  => $self->o('out'),
                  'extension' => '_short.csv' }
            ],
            -flow_into => {
                1 => { 'import_short_csv' => { 'csv_file' => '#ncr_file#' } },
            },
            -wait_for => [ 'import_long_csv' ]
        },

        {   -logic_name    => 'import_short_csv',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::ImportCsv',
        },

    ];
}

1;

