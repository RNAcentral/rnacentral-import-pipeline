#!/usr/bin/env perl

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
