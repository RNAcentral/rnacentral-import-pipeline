#!/usr/bin/env perl



=pod

=head1 NAME



=head1 SYNOPSIS

    Non-parallelized version of the data import pipeline.

    perl modules/RNAcentral/Import.pm --dir data --user=$user --password=$password --host=$host --sid=$sid --port=$port --out /path/to/temp/folder

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
    'host=s'     => \$opt->{'host'},
    'out=s'      => \$opt->{'out'}
);
&pod2usage if ( $location eq ''              or
                !defined($opt->{'user'})     or
                !defined($opt->{'password'}) or
                !defined($opt->{'sid'})      or
                !defined($opt->{'port'})     or
                !defined($opt->{'host'})     or
                !defined($opt->{'out'}));


# my $a = Bio::RNAcentral::InputFiles->new($opt);
# $a->process_folder($location);

test_merge($opt, $location, 'ncr');

exit;

my $b = Bio::RNAcentral::SqlldrImport->new($opt);
my @csvfiles = $a->list_folder($opt->{'out'}, 'csv');
for my $csvfile (@csvfiles) {
    $b->load_seq($csvfile);
}



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