#!/usr/bin/env perl

use strict;
use warnings;

use Cwd            ();
use File::Basename ();
BEGIN {
    my $root_dir = File::Basename::dirname( File::Basename::dirname( Cwd::realpath($0) ) );
    unshift @INC, "$root_dir/modules";
}

use Bio::RNAcentral::InputFiles;


my $location = '';
my $opt = {};

$opt->{'out'} = '/homes/apetrov/rnac-loader';

my $a = Bio::RNAcentral::InputFiles->new($opt);
my @files = $a->embl2csv('/homes/apetrov/rnac-loader/data/cum_std_mam_01_r117.ncr');

print join(',', @files), "\n";

