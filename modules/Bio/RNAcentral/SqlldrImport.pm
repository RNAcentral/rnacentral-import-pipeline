#!/usr/bin/env perl

package Bio::RNAcentral::SqlldrImport;

=pod

=head1 NAME

test

=head1 SYNOPSIS

=head1 DESCRIPTION

    Short sequences can be loaded very fast using direct path loading.
    Long sequences are stored as CLOBs in the Oracle database
    and cannot be loaded using direct path sqlldr strategy.
    As a result, long and short sequences are loaded separately.

    Creates temporary output files (sqlldr control files and sqlldr log files)
    in the directory specified by 

=head1 CONTACT


=cut

use strict;
use warnings;

use base ('Bio::RNAcentral::Base');


sub new {
    my ($class, $opt) = @_;

    # run parent constructor
    my $self = $class->SUPER::new();

    $self->{'job_id'}   = '';
    $self->{'seq_type'} = '';    

    $self->{'user'}     = $opt->{'user'};
    $self->{'password'} = $opt->{'password'};
    $self->{'sid'}      = $opt->{'sid'};
    $self->{'port'}     = $opt->{'port'};
    $self->{'host'}     = $opt->{'host'};

    # location of temporary output files
    if ( defined($opt->{'out'}) ) {
        $self->{'opt'}{'temp_dir'} = $opt->{'out'};
    }

    return $self;
}


=head2 load_seq

    Input: location of csv file with sequence data.
    Create control files for sqlldr, launch sqlldr, report mistakes.

=cut

sub load_seq {
    (my $self, my $file) = @_;

    # get file name without extension
    $file = File::Basename::fileparse($file, qr/\.[^.]*/);

    if ( $file =~ /(.*)_(long|short)/ ) {
        $self->{'job_id'}   = $1;
        $self->{'seq_type'} = $2;
    }    

    # create sqlldr control file
    $self->_make_ctl_file();

    $self->_delete_old_log_files();

    # launch sqlldr
    my $cmd = $self->_get_sqlldr_command();
    # $self->_run_sqlldr($cmd);

    # $self->{'logger'}->info($cmd);

    $self->_check_sqlldr_status();
}


=head2 _get_sqlldr_command

    Construct the sqlldr system command.

=cut

sub _get_sqlldr_command {
    my $self = shift;

    my $badfile = $self->get_output_filename($self->{'job_id'}, $self->{'seq_type'}, 'bad');
    my $logfile = $self->get_output_filename($self->{'job_id'}, $self->{'seq_type'}, 'log');
    my $ctlfile = $self->get_output_filename($self->{'job_id'}, $self->{'seq_type'}, 'ctl');

    my $command = 'sqlldr ' .
                   $self->{'user'} . '/' . $self->{'password'} . '@' .
                  '\"\(DESCRIPTION=\(ADDRESS=\(PROTOCOL=TCP\)' .
                  '\(HOST=' . $self->{'host'} . '\)' .
                  '\(PORT=' . $self->{'port'} . '\)\)' .
                  '\(CONNECT_DATA\=\(SERVICE_NAME=' . $self->{'sid'} . '\)\)\)\" ' .
                  'control=' . $ctlfile . ' ' .
                  'bad='     . $badfile . ' ' .
                  'log='     . $logfile;

    # conventional loading for lob files, direct loading for short sequences
    if ( $self->{'seq_type'} eq 'short' ) {
        $command .= ' direct=true parallel=true';
    }

    return $command;
}


=head2 _delete_old_log_files

    Delete any .bad and .log files left from previous runs.
    If they are not deleted, it may look as if the current run has an error too.

=cut

sub _delete_old_log_files {
    my $self = shift;

    my $badfile = $self->get_output_filename($self->{'job_id'}, $self->{'seq_type'}, 'bad');
    my $logfile = $self->get_output_filename($self->{'job_id'}, $self->{'seq_type'}, 'log');

    my @to_delete = ($badfile, $logfile);
    foreach my $file ( @to_delete ) {
        if ( -e $file ) {
            unlink $file or $self->{'logger'}->logwarn("Could not unlink $file: $!");
        }
    }
}

=head2

    Launch sqlldr and make sure that it runs successfully, log any errors.

=cut

sub _run_sqlldr {
    (my $self, my $command) = @_;
    $self->{'logger'}->info('Launching sqlldr');
    system($command) == 0
        or $self->{'logger'}->warn("Couldn't launch sqlldr\n. Command: $command\n Error: $!\n");
}


=head2

    Create a control file used by sqlldr.

=cut

sub _make_ctl_file {
    my $self = shift;

    my $ctlfile  = $self->get_output_filename($self->{'job_id'}, $self->{'seq_type'}, 'ctl');
    my $datafile = $self->get_output_filename($self->{'job_id'}, $self->{'seq_type'}, 'csv');

    open my $fh, '>', $ctlfile or die $!;

    print $fh <<CTL;
LOAD DATA
INFILE '$datafile' "str '\\n'"
APPEND
INTO TABLE $self->{'opt'}{'staging_table'}
FIELDS TERMINATED BY ','
(
  CRC64 char,
  LEN integer external,
CTL

    # different handling of long and short sequences
    if ( $self->{'seq_type'} eq 'short' ) {
        print $fh "  SEQ_SHORT char(" . $self->{'opt'}{'MAXSEQSHORT'} . "),\n";
    } elsif ( $self->{'seq_type'} eq 'long' ) {
        print $fh "  SEQ_LONG char("  . $self->{'opt'}{'MAXSEQLONG'}  . "),\n";
    } else {
        $self->{'logger'}->logdie("Wrong sequence sequence type");
    }

    print $fh <<CTL;
  AC char,
  VERSION integer external,
  TAXID integer external,
  MD5 char
)
CTL
    close $fh;
}


=head2 _check_sqlldr_status

    Entries rejected by the database are stored in the .bad file.
    Warn if bad file exists.

=cut

sub _check_sqlldr_status {
    my $self = shift;

    my $badfile = $self->get_output_filename($self->{'job_id'}, $self->{'seq_type'}, 'bad');

    if (-e $badfile) {
        $self->{'logger'}->warn("sqlldr import had errors");
    } else {
        $self->{'logger'}->info("No sqlldr errors");
    }

    # TODO: find error messages in the log file

}

1;