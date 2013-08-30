#!/usr/bin/env perl

=pod

=head1 NAME

    rnac_loader.pl

=head1 USAGE

    . ./config/hive_params
    perl scripts/rnac_loader.pl -in /path/to/ncr_files -out /path/to/temp_folder -user=$ORACLE_USER -password=$ORACLE_PASSWORD -host=$ORACLE_HOST -port=$ORACLE_PORT -sid=$ORACLE_SID

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
use Bio::RNAcentral::SqlldrImport;
use Bio::RNAcentral::OracleUpdate;
use Bio::RNAcentral::Embl2csv;
use Bio::RNAcentral::SqlldrImportReferences;


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
    'host=s'     => \$opt->{'host'}
);

&pod2usage if ( $location eq ''              or
                !defined($opt->{'user'})     or
                !defined($opt->{'password'}) or
                !defined($opt->{'sid'})      or
                !defined($opt->{'port'})     or
                !defined($opt->{'host'})     or
                !defined($opt->{'output_folder'}));

my $a = Bio::RNAcentral::InputFiles->new($opt);
my $e = Bio::RNAcentral::Embl2csv->new($opt);
my $b = Bio::RNAcentral::SqlldrImport->new($opt);
my $c = Bio::RNAcentral::OracleUpdate->new($opt);
my $d = Bio::RNAcentral::SqlldrImportReferences->new($opt);

# prepare staging table
$c->db_oracle_connect();
$c->truncate_staging_table();

# create output files in csv format
my @ncrfiles = $a->list_folder_recursive($location, 'ncr');
for my $ncrfile (@ncrfiles) {
    $e->embl2csv($ncrfile);
}

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

$d->load_all_references();

# launch plsql update
$c->run_pl_sql_update();
$c->db_oracle_disconnect();
