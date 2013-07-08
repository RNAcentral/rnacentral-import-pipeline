#!/usr/bin/env perl

# package RNAcentral::Import;
# TODO: complete package structure once it becomes clear how to integrate this with eHive

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

# lib folder contains SWISS:CRC64 and the modified version of BIO::SeqIO::embl.pm
# add it to the beginning of @INC
use FindBin qw($Bin);
use lib "$Bin/../../lib";

use Getopt::Long;
use Pod::Usage;
use Log::Log4perl;
use Bio::SeqIO;   # BioPerl is used for reading embl files
use SWISS::CRC64; # cyclic redundancy check
use Digest::MD5 qw(md5_hex);
use DBI;
use DBD::Oracle;


# Global parameters
our $OUTPUT        = '/homes/apetrov/rnac-loader/temp'; # TODO replace with a CLI option
our $STAGING_TABLE = 'load_rnacentral';
our $RELEASE_TABLE = 'rnc_release';
our $MAXSEQLONG    = 1000000;  # maximum length for long sequences stored as clobs
our $MAXSEQSHORT   = 4000;     # maximum length for short sequences stored as chars
our $EXTENSION     = 'ncr';    # look for .ncr files
our $SIZECUTOFF    = 10**7;    # 10 Mb, file size cutoff


main();

sub main {

    my $location = '';
    my @files = ();

    my $self = {};

    &GetOptions (
        'dir=s'      => \$location,
        'user=s'     => \$self->{'user'},
        'password=s' => \$self->{'password'},
        'sid=s'      => \$self->{'sid'},
        'port=i'     => \$self->{'port'},
        'host=s'     => \$self->{'host'}
    );
    &pod2usage if ( $location eq ''               or
                    !defined($self->{'user'})     or
                    !defined($self->{'password'}) or
                    !defined($self->{'sid'})      or
                    !defined($self->{'port'})     or
                    !defined($self->{'host'}) );

    $self->{'job_id'} = 'test_output'; # TODO: make dynamic, replace with the current filename
    $self->{'logger'} = initialize_logger();
    $self->{'dbh'}    = db_oracle_connect($self);

    # prepare the database
    truncate_staging_table($self);
    create_new_release($self);

    # get a list of all files
    my @original_files = list_files($self, $location, $EXTENSION);
    # my @all_files = ();
    # # split large files into small chunks and analyze them instead
    # foreach my $file (@original_files) {
    #     @all_files = (@all_files, file2chunks($self, $file));
    # }


    my @files = @original_files;

    # create csv files based on embl input files
    embl2csv(\@files, $self->{'job_id'});

    # load the data into the staging table using sqlldr
    load_seq($self, 'long');
    load_seq($self, 'short');

    # launch PL/SQL update
    run_pl_sql_update($self);

    db_oracle_disconnect($self);
}


=head2 initialize_logger

    Initialize the logger object and return it so it can be reused.
    TODO: switch to logging to file or both to screen and file

=cut

sub initialize_logger {

    # log4perl.appender.LOG1.filename  = mylog.log
    # log4perl.appender.LOG1.mode      = append

    my $log_conf = q(
       log4perl.rootLogger              = DEBUG, LOG1
       log4perl.appender.LOG1           = Log::Log4perl::Appender::Screen
       log4perl.appender.Screen.stderr  = 0
       log4perl.appender.LOG1.layout    = Log::Log4perl::Layout::PatternLayout
       log4perl.appender.LOG1.layout.ConversionPattern = %d %p %m %n
    );

    Log::Log4perl::init_once(\$log_conf);

    my $logger = Log::Log4perl->get_logger();

    $logger->info("Logger initialized");

    return $logger;
}


=head2 db_oracle_connect

    Establish a database connection and keep it open.

=cut

sub db_oracle_connect {
    my $self = shift;

    my $dsn = "dbi:Oracle:host=$self->{'host'};sid=$self->{'sid'};port=$self->{'port'}";
    my $dbh = DBI->connect($dsn, $self->{'user'}, $self->{'password'})
              or $self->{'logger'}->logdie( $DBI::errstr . "\n" );

    $self->{'logger'}->info("Connected to the database");
    return $dbh;
}

=head2 db_oracle_disconnect

    Close database connection.

=cut

sub db_oracle_disconnect {
    my $self = shift;

    if ( exists $self->{'dbh'} ) {
        $self->{'dbh'}->disconnect();
        $self->{'logger'}->info("Disconnected from the database");
    }
}


=head2 truncate_staging_table

    The staging table should be empty at the beginning of the import.

=cut

sub truncate_staging_table {
    my $self = shift;

    my $command = "TRUNCATE TABLE $STAGING_TABLE";
    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("Couldn't truncate the staging table");

    $self->{'logger'}->info("Staging table truncated");
}


=head2 create_new_release

    TODO: figure out how many staging tables are going to be used.
    TODO: replace placeholder parameters with real dates etc.

=cut

sub create_new_release {
    my $self = shift;

    my $command = <<SQL;
INSERT INTO $RELEASE_TABLE
VALUES(
  (SELECT count(*)+1 FROM rnc_release),
  1,
  '12-JUL-13',
  'F',
  'L',
  '12-JUL-13',
  'perl',
  'test',
  'N'
)
SQL

    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("Couldn't create new release");

    $self->{'logger'}->info("New release created");
}


=head2 load_seq

    Short sequences can be loaded very fast using direct path loading.
    Long sequences are stored as CLOBs in the Oracle database
    and cannot be loaded using direct path sqlldr strategy.
    As a result, long and short sequences are loaded separately.

=cut

sub load_seq {
    (my $self, my $seq_type) = @_;

    # create sqlldr control file
    _make_ctl_file($self, $seq_type);

    _delete_old_log_files($self, $seq_type);

    # launch sqlldr
    my $command = _get_sqlldr_command($self, $seq_type);
    _run_sqlldr($self, $command);

    _check_sqlldr_status($self, $seq_type);
}


=head2 _get_sqlldr_command

    Construct the sqlldr system command.

=cut

sub _get_sqlldr_command {
    (my $self, my $seq_type) = @_;

    my $badfile = get_filename($self->{'job_id'}, $seq_type, 'bad');
    my $logfile = get_filename($self->{'job_id'}, $seq_type, 'log');
    my $ctlfile = get_filename($self->{'job_id'}, $seq_type, 'ctl');

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
    if ( $seq_type eq 'short' ) {
        $command .= ' direct=true parallel=true';
    }

    return $command;
}


=head2 _delete_old_log_files

    Delete any .bad and .log files left from previous runs.
    If they are not deleted, it may look as if the current run has an error too.

=cut

sub _delete_old_log_files {
    (my $self, my $seq_type) = @_;

    my $badfile = get_filename($self->{'job_id'}, $seq_type, 'bad');
    my $logfile = get_filename($self->{'job_id'}, $seq_type, 'log');

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

    Generic function for constructing filenames for different types of files,
    such as ctl, csv, bad, and log files.

=cut

sub get_filename {

    (my $job_id, my $prefix, my $extension) = @_;

    return $OUTPUT . '/' . $job_id . '_' . $prefix . '.' . $extension;
}


=head2

    Create a control file used by sqlldr.

=cut

sub _make_ctl_file {
    (my $self, my $seq_type) = @_;

    my $ctlfile  = get_filename($self->{'job_id'}, $seq_type, 'ctl');
    my $datafile = get_filename($self->{'job_id'}, $seq_type, 'csv');

    my $fh = IO::File->new("> $ctlfile")
             or $self->{'logger'}->logdie("Couldn't open file $ctlfile");

    print $fh <<CTL;
LOAD DATA
INFILE '$datafile' "str '\\n'"
APPEND
INTO TABLE load_rnacentral
FIELDS TERMINATED BY ','
(
  CRC64 char,
  LEN integer external,
CTL

    # different handling of long and short sequences
    if ( $seq_type eq 'short' ) {
        print $fh "  SEQ_SHORT char($MAXSEQSHORT),\n";
    } elsif ( $seq_type eq 'long' ) {
        print $fh "  SEQ_LONG char($MAXSEQLONG),\n";
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
    $fh->close;
}


=head2 _check_sqlldr_status

    Entries rejected by the database are stored in the .bad file.
    Warn if bad file exists.

=cut

sub _check_sqlldr_status {

    my $self = shift @_;
    my $seq_type = shift @_;

    my $badfile = get_filename($self->{'job_id'}, $seq_type, 'bad');

    if (-e $badfile) {
        $self->{'logger'}->warn("sqlldr import had errors");
    } else {
        $self->{'logger'}->info("No sqlldr errors");
    }

    # TODO: find error messages in the log file

}


=head2

    Read a list of embl files and create two csv files,
    one with short, another with long sequences.

=cut

sub embl2csv {

    my @files  = @{$_[0]};
    my $job_id = $_[1];

    my $fname_long  = get_filename($job_id, 'long', 'csv');
    my $fname_short = get_filename($job_id, 'short', 'csv');

    # open output file
    my $fh_long  = IO::File->new("> $fname_long")  or die("Couldn't open file $fname_long");
    my $fh_short = IO::File->new("> $fname_short") or die("Couldn't open file $fname_short");

    # loop over files
    foreach my $file (@files) {

        my $stream = Bio::SeqIO->new(-file => $file, -format => 'EMBL');

        # reset counters
        my $i = 0;
        my $seqs_long  = 0;
        my $seqs_short = 0;

        # loop over records
        while ( (my $seq = $stream->next_seq()) ) {

            $i++;
            # print "$i\n";

            # reset all values to avoid accidental reuse of data from previous records
            my $ac      = '';
            my $length  = 0;
            my $version = 0;
            my $md5     = '';
            my $crc64   = '';
            my $taxid   = '';
            my $isLong  = -1;
            my $data    = '';
            my $sequence = '';

            # get new data
            $ac       = $seq->display_id;
            $length   = $seq->length;
            $version  = $seq->seq_version;
            $sequence = $seq->seq;

            if ($length > $MAXSEQSHORT) {
                $isLong = 1;
                $seqs_long++;
            } else {
                $seqs_short++;
                $isLong = 0;
            }

            $md5   = md5_hex($seq->seq);
            $crc64 = SWISS::CRC64::crc64($seq->seq);

            for my $feat_object ($seq->get_SeqFeatures) {
                if ( $feat_object->primary_tag eq 'source' ) {
                    for my $tag ($feat_object->get_all_tags) {
                        for my $value ($feat_object->get_tag_values($tag)) {
                            if ($tag eq 'db_xref') {
                                (my $text, $taxid) = split(':', $value);
                                last;
                            }
                        }
                    }
                    last; # found taxid, stop looking
                }
            }

            $data = join(',', ($crc64, $length, $sequence, $ac, $version, $taxid, $md5)) . "\n";

            # print the data in different files depending on sequence length
            if ( $isLong == 1 ) {
                print $fh_long $data;
            } elsif ( $isLong == 0 ) {
                print $fh_short $data;
            } else {
                print 'Error';
            }
        }
    }

    $fh_short->close;
    $fh_long->close;
}


=head2 list_files

    Get all files with the specified extension from the input directory

=cut

sub list_files {

    (my $self, my $dir, my $extension) = @_;

    opendir(DIR, $dir) or die $!;

    my @files = ();
    while (my $file = readdir(DIR)) {
        next if ($file =~ m/^\./);
        next if ($file !~ m/\.$extension$/);

        push @files, $dir . '/' . $file;
    }

    closedir(DIR);

    $self->{'logger'}->info("Found " . scalar(@files) . " .$extension files");

    return @files;
}


=head2 file2chunks

    Some input files are very large, so it's faster to split them into chunks for processing.
    Split large files into chunks and return an array of filenames.
    Return an array with the original file if it's small.

=cut

sub file2chunks {

    (my $self, my $filename) = @_;

    my $size = -s $filename;
    if ( $size < $SIZECUTOFF ) {
        return ($filename);
    }

    $self->{'logger'}->info("File $filename is $size bytes, will split it into chunks");

    local $/ = "//\n";

    my $i      = 1;
    my $text   = '';
    my @chunks = ();
    my $chunk  = '';
    (my $path, my $extension) = split('\.', $filename);

    use bytes; # to calculate length in bytes

    open (INFILE, $filename);
    while (<INFILE>) {
        $text .= $_;
        if ( length($text) > $SIZECUTOFF ) {
            $chunk = $path . '_chunk' . $i . '.' . $extension;
            push @chunks, $chunk;
            _print_to_file($chunk, $text);
            $i++;
            $text = '';
            $self->{'logger'}->info("Created file $chunk");
        }
    }
    close INFILE;

    # left over sequences
    if ( length($text) > 0 ) {
        $chunk = $path . '_chunk' . $i . '.' . $extension;
        _print_to_file($chunk, $text);
        $self->{'logger'}->info("Created file $chunk");
    }

    return @chunks;
}


sub _print_to_file {
    open (OUTFILE, "> $_[0]") or die("Couldn't open file");
    print OUTFILE $_[1];
    close OUTFILE;
}


=head2 run_pl_sql_update

    Update Oracle database, import sequences, create new xrefs.

=cut

sub run_pl_sql_update {
    my $self = shift;

    my $command = <<PLSQL;
BEGIN
  RNC_UPDATE.NEW_UPDATE();
END;
PLSQL

    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("PL/SQL update failed");

    $self->{'logger'}->info("PL/SQL update complete");
}
