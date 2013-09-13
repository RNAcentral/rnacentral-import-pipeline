#!/usr/bin/env perl

package Bio::RNAcentral::SqlldrImportCompositeIds;

=pod

=head1 NAME

    Bio::RNAcentral::SqlldrImportCompositeIds

=head1 DESCRIPTION



=cut

use strict;
use warnings;

use File::Spec;

use base ('Bio::RNAcentral::SqlldrImportBase');


sub new {
    my ($class, $opt, $path, $prefix) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt, $path, $prefix);

    return $self;
}


=head2 load_all_references

    Main subroutine for loading the references.

=cut

sub load_all {
    my $self = shift;

    $self->{'logger'}->info("Loading composite ids");

    # create one sqlldr control file
    $self->_make_ctl_file();

    # concatenate all csv files
    my $cmd = "cat $self->{'local'}{'path'}/*.csv > $self->{'local'}{'allfile'}";
    $self->{'logger'}->info("Creating a single datafile with composite ids: $cmd");
    system($cmd);

    # clean up old files if present
    $self->_delete_old_log_files();

    # prepare sqlldr command
    $cmd = $self->_get_sqlldr_command();

    # run sqlldr
    $self->{'logger'}->info("Importing $self->{'local'}{'allfile'}");
    my $problems = $self->_run_sqlldr($cmd);

    # clean up if no errors and no problems in sqlldr
    unless ( $self->_errors_found() or $problems ) {
        unlink $self->{'local'}{'logfile'}, $self->{'local'}{'badfile'};
    }
}


=head2 _make_ctl_file

    Create a control file used by sqlldr.

=cut

sub _make_ctl_file {
    my $self = shift;

    open my $fh, '>', $self->{'local'}{'ctlfile'} or die $!;

    print $fh <<CTL;
LOAD DATA
INFILE "str '\\n'"
APPEND
INTO TABLE $self->{'opt'}{'comp_id_table'}
FIELDS TERMINATED BY ',' enclosed by '"'
(
    COMPOSITE_ID char,
    AC char,
    DATABASE char,
    OPTIONAL_ID char,
    EXTERNAL_ID char
)
CTL
    close $fh;
}


=head2 update

    Import the data and perform pl/sql update.

=cut

sub update {
    my $self = shift;

    $self->db_oracle_connect();
    $self->truncate_table($self->{'opt'}{'comp_id_table'});
    $self->load_all();
    $self->update_composite_ids();
    $self->db_oracle_disconnect();
}


1;