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


my $location = '';
my $opt = {};

# all options are mandatory
&GetOptions (
    'in=s'       => \$location,          # input data location
    'out=s'      => \$opt->{'out'},      # location of output temp files
    'user=s'     => \$opt->{'user'},     # Oracle connection details
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
                !defined($opt->{'out'}));

my $a = Bio::RNAcentral::InputFiles->new($opt);
my $b = Bio::RNAcentral::SqlldrImport->new($opt);
my $c = Bio::RNAcentral::OracleUpdate->new($opt);

# prepare staging table
$c->truncate_staging_table();


my @folders = ($location); # list of folders with ncr files

for my $folder (@folders) {
    # prepare csv files for sqlldr
    $a->process_folder($folder);

    # load data into the staging table
    my @csvfiles = $a->list_folder($opt->{'out'}, 'csv');
    for my $csvfile (@csvfiles) {
        $b->load_seq($csvfile);
    }
}


# launch plsql update
# $c->create_new_release();
# $c->run_pl_sql_update();
$c->db_oracle_disconnect();



sub test_merge {

    my ($opt, $location, $extension) = @_;
    my $self = Bio::RNAcentral::InputFiles->new($opt);
    my ($size_ref, $ordered_files_ref) = $self->list_folder_order_by_size($location, $extension);

    my $groups_ref = $self->group_files($size_ref, $ordered_files_ref);

    use Data::Dumper;
    print Dumper($groups_ref);

    # verify that all files were found
    my @arr = @{$groups_ref};
    my $total = 0;
    foreach my $group (@arr) {
        if ($group) {
            $total += scalar @$group;
        }
    }

    print $total, "\n";

}