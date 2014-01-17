#!/usr/bin/env perl

=pod

=head1 NAME

    update_bibliographic_data.pl

=head1 USAGE

    . ./config/hive_params
    perl scripts/update_bibliographic_data.pl -user=$ORACLE_USER -password=$ORACLE_PASSWORD -host=$ORACLE_HOST -port=$ORACLE_PORT -sid=$ORACLE_SID

=head1 DESCRIPTION

    Import missing bibliographic data.

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

use Bio::RNAcentral::UpdateRefsUsingEuropePMC;


my $opt = {};

# all options are mandatory
&GetOptions (
    'user=s'     => \$opt->{'user'},
    'password=s' => \$opt->{'password'},
    'sid=s'      => \$opt->{'sid'},
    'port=i'     => \$opt->{'port'},
    'host=s'     => \$opt->{'host'},
);

&pod2usage if ( !defined($opt->{'user'})     or
                !defined($opt->{'password'}) or
                !defined($opt->{'sid'})      or
                !defined($opt->{'port'})     or
                !defined($opt->{'host'}));

$a = Bio::RNAcentral::UpdateRefsUsingEuropePMC->new($opt);
$a->db_oracle_connect();
$a->update_bibliographic_info();
$a->db_oracle_disconnect();
