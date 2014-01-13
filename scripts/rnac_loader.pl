#!/usr/bin/env perl
# Copyright [2009-2014] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


=pod

=head1 NAME

    rnac_loader.pl

=head1 USAGE

    . ./config/hive_params
    perl scripts/rnac_loader.pl -in /path/to/ncr_files -out /path/to/temp_folder -user=$ORACLE_USER -password=$ORACLE_PASSWORD -host=$ORACLE_HOST -port=$ORACLE_PORT -sid=$ORACLE_SID -release_type=F/I

=head1 DESCRIPTION

    Non-parallelized version of the data import pipeline.

=cut

use strict;
use warnings;

use Cwd            ();
use File::Basename ();
BEGIN {
    my $root_dir = File::Basename::dirname( File::Basename::dirname( Cwd::realpath($0) ) );
    unshift @INC, "$root_dir/modules";
}

use Getopt::Long;
use Pod::Usage;

use Bio::RNAcentral::InputFiles;
use Bio::RNAcentral::OracleUpdate;
use Bio::RNAcentral::Embl2csv;
use Bio::RNAcentral::SqlldrImportSequences;
use Bio::RNAcentral::SqlldrImportReferences;
use Bio::RNAcentral::SqlldrImportAccessionInfo;
use Bio::RNAcentral::SqlldrImportCompositeIds;


my $location = '';
my $opt = {};

# all options are mandatory
&GetOptions (
    'in=s'       => \$location,               # input data location
    'out=s'      => \$opt->{'output_folder'}, # location of output temp files
    'user=s'     => \$opt->{'user'},          # Oracle connection details
    'password=s' => \$opt->{'password'},
    'sid=s'      => \$opt->{'sid'},
    'port=i'     => \$opt->{'port'},
    'host=s'     => \$opt->{'host'},
    'release_type=s' => \$opt->{'release_type'}, # release type (full or incremental)
);

&pod2usage if ( $location eq ''              or
                !defined($opt->{'user'})     or
                !defined($opt->{'password'}) or
                !defined($opt->{'sid'})      or
                !defined($opt->{'port'})     or
                !defined($opt->{'host'})     or
                !defined($opt->{'output_folder'}));

# Full release is default
if ( !defined($opt->{'release_type'}) ) {
    $opt->{'release_type'} = 'F';
}

my $a = Bio::RNAcentral::InputFiles->new($opt);
my $e = Bio::RNAcentral::Embl2csv->new($opt);
my $b = Bio::RNAcentral::SqlldrImportSequences->new($opt);
my $c = Bio::RNAcentral::OracleUpdate->new($opt);

# prepare staging table
$c->db_oracle_connect();
$c->truncate_table($c->{'opt'}{'staging_table'});

# create output files in csv format
my @ncrfiles = $a->list_folder_recursive($location, 'ncr');
for my $ncrfile (@ncrfiles) {
    $e->embl2csv($ncrfile);
}

# load information about non-coding accessions
my $f = Bio::RNAcentral::SqlldrImportAccessionInfo->new($opt, 'ac_info');
$f->update();

# load information about composite non-coding ids
my $g = Bio::RNAcentral::SqlldrImportCompositeIds->new($opt, 'composite_ids');
$g->update();

# load literature references
my $d = Bio::RNAcentral::SqlldrImportReferences->new($opt, 'refs');
$d->update();

# load long sequences
$b->make_ctl_files();
my @csvfiles = $a->list_folder($a->get_long_folder_path(), 'csv');
for my $csvfile (@csvfiles) {
    $b->load_seq($csvfile);
}

# load short sequences
@csvfiles = $a->list_folder($a->get_short_folder_path(), 'csv');
for my $csvfile (@csvfiles) {
    $b->load_seq($csvfile);
}

# launch plsql update
$c->run_pl_sql_update();
$c->db_oracle_disconnect();
