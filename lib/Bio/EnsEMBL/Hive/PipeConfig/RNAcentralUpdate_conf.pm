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

    Bio::EnsEMBL::Hive::PipeConfig::RNAcentralUpdate_conf

=head1 DESCRIPTION

    Hive pipeline configuration file for the RNAcentral data import pipeline.

=cut


package Bio::EnsEMBL::Hive::PipeConfig::RNAcentralUpdate_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Use the command line parameters to configure the pipeline.

=cut

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },    # inherit other stuff from the base class

        'pipeline_name' => 'rnacentral_staging', # name used by the beekeeper to prefix job names on the farm
        'hive_force_init' => 1,                  # always recreate the hive database

        # get the command line options
        'pipeline_db'   => {
            -host   => $self->o('host'),
            -port   => $self->o('port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -dbname => $self->o('pipeline_name'),
        },

    };
}


=head2 pipeline_create_commands

    Hive internal processing.

=cut

sub pipeline_create_commands {
    my ($self) = @_;
    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation

    ];
}


=head2 pipeline_wide_parameters

    Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
    The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
    Please see existing PipeConfig modules for examples.

=cut

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class

        # store command line parameters
        'output_folder'   => $self->o('output_folder'),
        'oracle-user'     => $self->o('oracle-user'),
        'oracle-password' => $self->o('oracle-password'),
        'oracle-sid'      => $self->o('oracle-sid'),
        'oracle-port'     => $self->o('oracle-port'),
        'oracle-host'     => $self->o('oracle-host'),
        'release_type'    => $self->o('release_type'), # incremental or full
    };
}


=head2 pipeline_analyses

    Main logic of the pipeline.

=cut

sub pipeline_analyses {
    my ($self) = @_;

    return [
        # truncate staging table
        {   -logic_name => 'truncate_staging_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::TruncateStagingTable',
            -analysis_capacity  => 1,
            -input_ids  => [
                { 'id' => 1 }
            ],
        },

        # # get ncr files from a local folder
        # {   -logic_name => 'get_ncr_files',
        #     -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetFiles',
        #     -analysis_capacity  => 1,
        #     -input_ids  => [
        #         { 'location'  => '',
        #           'extension' => 'ncr' }
        #     ],
        #     -flow_into => {
        #         1 => { 'create_csv_files' => { 'ncr_file' => '#ncr_file#' } },
        #     },
        #     -wait_for => [ 'truncate_staging_table' ]
        # },

        # list all ncr files from the ftp site
        {   -logic_name => 'get_ncr_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetNcProduct',
            -input_ids  => [
                { 'release_type'  => $self->o('release_type') }
            ],
            -analysis_capacity  =>  1,
            -flow_into => {
                1 => { 'create_csv_files' => { 'ncr_file' => '#ncr_file#' } },
            },
        },

        # prepare csv files for sql loader
        {   -logic_name    => 'create_csv_files',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::CreateCsvFiles',
        },

        # load literature references
        {   -logic_name => 'load_references',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::LoadReferences',
            -analysis_capacity  => 1,
            -input_ids  => [
                { 'id' => 1 }
            ],
            -wait_for => [ 'get_ncr_files', 'create_csv_files' ]
        },

        # load accession info
        {   -logic_name => 'load_accession_info',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::LoadAccessionInfo',
            -analysis_capacity  => 1,
            -input_ids  => [
                { 'id' => 1 }
            ],
            -wait_for => [ 'get_ncr_files', 'create_csv_files' ]
        },


        # list all csv files with sequences longer than 4000 characters.
        {   -logic_name => 'get_long_csv_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetFiles',
            -analysis_capacity  => 1,
            -input_ids  => [
                { 'location'  => $self->o('output_folder') . '/long',
                  'extension' => 'csv' }
            ],
            -flow_into => {
                1 => { 'import_long_csv' => { 'csv_file' => '#ncr_file#' } },
            },
            -wait_for => [ 'get_ncr_files', 'create_csv_files' ]
        },

        # long sequences need to be imported sequentially, because they are stored as clobs
        # and parallel loading for clobs is not supported in Oracle.
        {   -logic_name    => 'import_long_csv',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::LoadSequences',
            -analysis_capacity  =>  1,
        },

        # list all csv files with sequences shorter than 4000 characters
        {   -logic_name => 'get_short_csv_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::GetFiles',
            -analysis_capacity  => 1,
            -input_ids  => [
                { 'location'  => $self->o('output_folder') . '/short',
                  'extension' => 'csv' }
            ],
            -flow_into => {
                1 => { 'import_short_csv' => { 'csv_file' => '#ncr_file#' } },
            },
            -wait_for => [ 'import_long_csv' ]
        },

        # import all csv files with short sequences
        # highly parallel step
        {   -logic_name    => 'import_short_csv',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::LoadSequences',
        },

        # launch PL/SQL update
        {   -logic_name => 'launch_plsql_update',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::RNAcentral::LaunchPlsqlUpdate',
            -analysis_capacity  => 1,
            -input_ids  => [
                { 'release_type' => $self->o('release_type') }
            ],
            -wait_for => [ 'import_short_csv' ],
            -max_retry_count => 0,
        },

    ];
}

1;

