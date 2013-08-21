#!/usr/bin/env perl

package Bio::RNAcentral::SqlldrImportReferences;

=pod

=head1 NAME

    Bio::RNAcentral::SqlldrImportReferences

=head1 DESCRIPTION

    Import literature references into the RNAcentral Oracle database.
    Since some fields (e.g. authors) are longer than 4000 characters, they are stored as clobs.
    Parallel load option is not allowed when loading lob columns (SQL*Loader-971),
    so the data are loaded sequentially.

=cut

use strict;
use warnings;

use File::Spec;

use base ('Bio::RNAcentral::Base');


sub new {
    my ($class, $opt) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt);

    my $refs_path = $self->get_refs_path();

    $self->{'refs'} = {
        refs_path => $refs_path,
        ctlfile   => File::Spec->catfile($refs_path, 'refs.ctl'),
        allrefs   => File::Spec->catfile($refs_path, 'allrefs.dat'),
        badfile   => File::Spec->catfile($refs_path, 'allrefs.bad'),
        logfile   => File::Spec->catfile($refs_path, 'allrefs.log'),
    };

    return $self;
}


=head2 load_all_references

    Main subroutine for loading the references.

=cut

sub load_all_references {
    my $self = shift;

    $self->{'logger'}->info("Loading references");

    # create one sqlldr control file
    $self->_make_ctl_file();

    # concatenate all csv files
    my $cmd = "cat $self->{'refs'}{'refs_path'}/*.csv > $self->{'refs'}{'allrefs'}";
    $self->{'logger'}->info("Creating a single datafile with references: $cmd");
    system($cmd);

    # clean up old files if present
    $self->_delete_old_log_files();

    # prepare sqlldr command
    $cmd = $self->_get_sqlldr_command($self->{'refs'}{'allrefs'});

    # run sqlldr
    $self->{'logger'}->info("Importing $self->{'refs'}{'allrefs'}");
    my $problems = $self->_run_sqlldr($cmd);

    # clean up if no errors and no problems in sqlldr
    unless ( $self->_errors_found() or $problems ) {
        unlink $self->{'refs'}{'allrefs'}, $self->{'refs'}{'logfile'};
    }
}


=head2 _get_sqlldr_command

    Construct the sqlldr system command.

=cut

sub _get_sqlldr_command {
    my $self = shift;

    my $command = 'sqlldr ' .
                   $self->{'user'} . '/' . $self->{'password'} . '@' .
                  '\"\(DESCRIPTION=\(ADDRESS=\(PROTOCOL=TCP\)' .
                  '\(HOST=' . $self->{'host'} . '\)' .
                  '\(PORT=' . $self->{'port'} . '\)\)' .
                  '\(CONNECT_DATA\=\(SERVICE_NAME=' . $self->{'sid'} . '\)\)\)\" ' .
                  'control=' . $self->{'refs'}{'ctlfile'} . ' ' .
                  'bad='     . $self->{'refs'}{'badfile'} . ' ' .
                  'log='     . $self->{'refs'}{'logfile'} . ' ' .
                  'direct=true errors=99999999 data=' . $self->{'refs'}{'allrefs'};
    return $command;
}


=head2 _delete_old_log_files

    Delete any .bad and .log files left from previous runs.
    If they are not deleted, it may look as if the current run had an error and not the old one.

=cut

sub _delete_old_log_files {
    my $self = shift;
    my @to_delete = ($self->{'refs'}{'badfile'}, $self->{'refs'}{'logfile'});
    unlink @to_delete;
}

=head2

    Launch sqlldr and make sure that it runs successfully, log any errors.

=cut

sub _run_sqlldr {
    (my $self, my $command) = @_;
    $self->{'logger'}->info('Launching sqlldr');

    my $status = system($command); # 0 on success
    unless ( $status == 0 ) {
        $self->{'logger'}->logwarn("Couldn't launch sqlldr\n. Command: $command\n Error: $!\n");
    }

    return $status;
}


=head2

    Create a control file used by sqlldr.

=cut

sub _make_ctl_file {
    my $self = shift;

    open my $fh, '>', $self->{'refs'}{'ctlfile'} or die $!;

    print $fh <<CTL;
LOAD DATA
INFILE "str '\\n'"
APPEND
INTO TABLE $self->{'opt'}{'references_table'}
FIELDS TERMINATED BY ',' enclosed by '"'
(
    MD5 char,
    AUTHORS_MD5 char,
    AUTHORS char(1000000),
    LOCATION char(4000),
    TITLE char(4000),
    PUBMED char,
    DOI char,
    PUBLISHER char,
    EDITORS char
)
CTL
    close $fh;
}


=head2 _errors_found

    Warn if bad file exists.
    Bad file contains the entries rejected by the database.

=cut

sub _errors_found {
    my $self = shift;

    # TODO: find error messages in the log file

    if (-e $self->{'refs'}{'badfile'}) {
        $self->{'logger'}->logwarn("sqlldr import had errors, check $self->{'refs'}{'badfile'}");
        return 1;
    } else {
        $self->{'logger'}->info("No bad file");
        return 0;
    }
}


1;