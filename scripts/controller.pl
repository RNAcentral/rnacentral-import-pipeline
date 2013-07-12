#!/usr/bin/env perl



=pod

=head1 NAME

test

=head1 SYNOPSIS

    perl modules/RNAcentral/Import.pm --dir data --user=$user --password=$password --host=$host --sid=$sid --port=$port

    data is the folder with input files in embl format.

=head1 DESCRIPTION

=head1 CONTACT


=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use Bio::RNAcentral::InputFiles;
use Bio::RNAcentral::SqlldrImport;


my $location = '';

my $opt = {};

&GetOptions (
    'dir=s'      => \$location,
    'user=s'     => \$opt->{'user'},
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
                !defined($opt->{'host'}) );


# my $a = Bio::RNAcentral::InputFiles->new();

# $a->process_folder($location);


my $b = Bio::RNAcentral::SqlldrImport->new($opt);
$b->load_seq('/Users/apetrov/Desktop/ensembl_main/rnac-loader/temp/rel_est_vrt_r115_long.csv');
