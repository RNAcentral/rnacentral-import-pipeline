#!/usr/bin/env perl

package Bio::RNAcentral::SqlldrImportReferences;

=pod

=head1 NAME

    Bio::RNAcentral::SqlldrImportReferences

=head1 DESCRIPTION


=cut

use strict;
use warnings;

use File::Spec;

use base ('Bio::RNAcentral::Base');


sub new {
    my ($class, $opt) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt);

    $self->{'refs'} = {
        refs_path => $self->get_refs_path(),
        ctlfile   => '',
        badfile   => '',
        logfile   => '',
    };

    # get ctl filename
    $self->{'refs'}{'ctlfile'} = File::Spec->catfile($self->{'output_folder'}, 'refs.ctl');

    return $self;
}


=head2 load_all_references




=cut

sub load_all_references {
    my $self = shift;

    $self->{'logger'}->info("Loading references");

    # create one sqlldr control file
    $self->_make_ctl_file();

    # get all csv files with references
    my $csv_files = File::Spec->catfile($self->{'refs'}{'refs_path'}, '*.csv');
    my @files = glob "$csv_files";
    $self->{'logger'}->info('Found ' . scalar(@files) . ' files with references');

    foreach my $file (@files) {
        $self->{'logger'}->info("Importing $file");

        # get file name without extension
        my $job_id = File::Basename::fileparse($file, qr/\.[^.]*/);

        # initialize all filenames for this job_id
        $self->_set_filenames($job_id);
        $self->_delete_old_log_files();

        # prepare sqlldr command
        my $cmd = $self->_get_sqlldr_command($file);

        # run sqlldr
        my $problems = $self->_run_sqlldr($cmd);

        # clean up if no errors and no problems in sqlldr
        unless ( $self->_errors_found() or $problems ) {
            # unlink $file, $self->{'refs'}{'logfile'};
        }
    }
}


=head2 _set_filenames

    Initialize all filenames once for easy access.

=cut

sub _set_filenames {
    (my $self, my $job_id) = @_;

    $self->{'refs'}{'logfile'} = File::Spec->catfile($self->{'refs'}{'refs_path'}, $job_id . '.log');
    $self->{'refs'}{'badfile'} = File::Spec->catfile($self->{'refs'}{'refs_path'}, $job_id . '.bad');
}


=head2 _get_sqlldr_command

    Construct the sqlldr system command.

=cut

sub _get_sqlldr_command {
    (my $self, my $file) = @_;

    my $command = 'sqlldr ' .
                   $self->{'user'} . '/' . $self->{'password'} . '@' .
                  '\"\(DESCRIPTION=\(ADDRESS=\(PROTOCOL=TCP\)' .
                  '\(HOST=' . $self->{'host'} . '\)' .
                  '\(PORT=' . $self->{'port'} . '\)\)' .
                  '\(CONNECT_DATA\=\(SERVICE_NAME=' . $self->{'sid'} . '\)\)\)\" ' .
                  'control=' . $self->{'refs'}{'ctlfile'} . ' ' .
                  'bad='     . $self->{'refs'}{'badfile'} . ' ' .
                  'log='     . $self->{'refs'}{'logfile'} . ' ' .
                  'direct=true data=' . $file;
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

    Entries rejected by the database are stored in the .bad file.
    Warn if bad file exists.
    Discard file is not specified in control files, so it's never created.

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