#!/usr/bin/env perl

# package RNAcentral::Import;

#

use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/../../lib"; # add lib folder to the beginning of @INC, this lets us run modified version of BIO::SeqIO::embl.pm

use Getopt::Long;
use Log::Log4perl;
use Bio::SeqIO;
use Digest::MD5 qw(md5_hex);
use SWISS::CRC64;
use DBI;
use DBD::Oracle;


our $OUTPUT = '/homes/apetrov/rnac-loader/temp'; # TODO replace with a CLI option
our $STAGING_TABLE = 'load_rnacentral';

main();

sub main {

    my $dir = '';
    my @files = ();

    my $self = {};
    $self->{'logger'} = initialize_logger();
    $self->{'job_id'} = 'test_output';

    &GetOptions (
        'dir=s'      => \$dir,
        'user=s'     => \$self->{'user'},
        'password=s' => \$self->{'password'},
        'sid=s'      => \$self->{'sid'},
        'port=i'     => \$self->{'port'},
        'host=s'     => \$self->{'host'}
    );
    &Usage if ($dir eq '');

    # $self->{'job_id'} = 'test_output';
    # $self->{'user'} = 'RNACEN';
    # $self->{'password'} ='necanr13';
    # $self->{'sid'}  ='RNCDEV';
    # $self->{'port'} =1541;
    # $self->{'host'} ='ora-vm-014.ebi.ac.uk';


    # prepare the database: truncate staging table, disable indexes
    prepare_staging_table($self);
    # create_new_release($self);

    # create csv files
    @files = get_embl_files($self, $dir);
    embl2csv(\@files, $self->{'job_id'});

    # load the data into the staging table
    load_seq($self, 'long');
    load_seq($self, 'short');

    # launch PL/SQL update


}


### SUBS ###


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


sub run_pl_sql_update {
    my $self = shift;

    my $dsn = "dbi:Oracle:host=$self->{'host'};sid=$self->{'sid'};port=$self->{'port'}";
    my $dbh = DBI->connect($dsn, $self->{'user'}, $self->{'password'})
              or $self->{'logger'}->logdie( $DBI::errstr . "\n" );

    my $command = <<PLSQL;
BEGIN
  RNC_UPDATE.NEW_UPDATE();
END;
PLSQL

    $dbh->do($command) or $self->{'logger'}->logdie("PL/SQL update failed");
    $dbh->disconnect();

    $self->{'logger'}->info("PL/SQL update complete");
}


sub prepare_staging_table {
    my $self = shift;

    my $dsn = "dbi:Oracle:host=$self->{'host'};sid=$self->{'sid'};port=$self->{'port'}";
    my $dbh = DBI->connect($dsn, $self->{'user'}, $self->{'password'})
              or $self->{'logger'}->logdie( $DBI::errstr . "\n" );

    my $command = "TRUNCATE TABLE $STAGING_TABLE";

    $dbh->do($command) or $self->{'logger'}->logdie("Couldn't truncate the staging table");
    $dbh->disconnect();

    $self->{'logger'}->info("Staging table ready");
}


# todo use a single dbh connection
sub create_new_release {
    my $self = shift;

    my $dsn = "dbi:Oracle:host=$self->{'host'};sid=$self->{'sid'};port=$self->{'port'}";
    my $dbh = DBI->connect($dsn, $self->{'user'}, $self->{'password'})
              or $self->{'logger'}->logdie( $DBI::errstr . "\n" );

    my $command = <<SQL;
INSERT INTO rnc_release
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

    $self->{'logger'}->info($command);

    $dbh->do($command) or $self->{'logger'}->logdie("Couldn't create new release");
    $dbh->disconnect();

    $self->{'logger'}->info("New release created");
}


sub load_seq {
    (my $self, my $seq_type) = @_;

    # create sqlldr control file
    make_ctl_file($self, $seq_type);

    # delete_old_log_files;

    # launch sqlldr
    my $command = get_sqlldr_command($self, $seq_type);
    system_call($command);

    check_sqlldr_status($self, $seq_type);
}


sub get_sqlldr_command {
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

    # direct conventional loading for lob files
    if ( $seq_type eq 'short' ) {
        $command .= ' direct=true parallel=true';
    }

    $self->{'logger'}->info($command);

    return $command;
}


# sub delete_old_log_files {
#     # delete bad and log files from previous runs
#     my @to_delete = ($badfile, $logfile);
#     foreach my $file ( @to_delete ) {
#         if ( -e $file ) {
#             unlink $file or warn "Could not unlink $file: $!";
#         }
#     }
# }


sub system_call {
    system(shift) == 0 or print STDERR "Couldn't launch sqlldr: $!\n";
}


sub get_filename {

    (my $job_id, my $prefix, my $extension) = @_;

    return $OUTPUT . '/' . $job_id . '_' . $prefix . '.' . $extension;
}


sub make_ctl_file {
    (my $self, my $seq_type) = @_;

    my $MAXSEQLONG = 1000000;

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

    if ( $seq_type eq 'short' ) {
        print $fh "  SEQ_SHORT char(4000),\n";
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


sub check_sqlldr_status {

    my $self = shift @_;
    my $seq_type = shift @_;

    my $badfile = get_filename($self->{'job_id'}, $seq_type, 'bad');

    # check if bad file exists
    if (-e $badfile) {
        print "Bad file exists, there were mistakes\n";
    } else {
        print "No mistakes\n";
    }

    # find error messages in the log file

}


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

            if ($length > 4000) {
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


sub get_embl_files {

    (my $self, my $dir) = @_;

    opendir(DIR, $dir) or die $!;

    my @files = ();
    while (my $file = readdir(DIR)) {
        next if ($file =~ m/^\./);
        next if ($file !~ m/\.ncr$/);

        push @files, $dir . '/' . $file;
    }

    closedir(DIR);

    $self->{'logger'}->info("Found " . scalar(@files) . " files\n");

    return @files;
}


sub Usage {
    print <<EOF
 $0 --dir <path to input files>

        Mandatory
            --dir     Location of input files
EOF
;
    exit(1);
}

__DATA__

== pod





== cut


