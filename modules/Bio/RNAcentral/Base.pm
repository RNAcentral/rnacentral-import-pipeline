=pod 

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION

=cut

package Bio::RNAcentral::Base;

use strict;

use Log::Log4perl;


our $OUTPUT = '/Users/apetrov/Desktop/ensembl_main/rnac-loader/temp'; # TODO replace with a CLI option


=head2 new

	Constructor, initializes logging.

=cut

sub new {
    my ($class, @args) = @_;

    my $self = bless {}, $class;

    $self->{'logger'} = initialize_logger();
    $self->{'opt'} = default_options();    

    return $self;
}


sub default_options {
    
    return {

        'file_extension'   => 'ncr',         # look for .ncr files
        'file_size_cutoff' => 50 * 10**6,    # 50 Mb, file size cutoff

        'staging_table' => 'load_rnacentral', # staging table
        'release_table' => 'rnc_release',     # keeps track of all RNAcentral releases
        'temp_dir'      => '/Users/apetrov/Desktop/ensembl_main/rnac-loader/temp2',       # default location of temporary output files, can be specified on command line
		'MAXSEQLONG'    => 1000000,  # maximum length for long sequences stored as clobs
		'MAXSEQSHORT'   => 4000,     # maximum length for short sequences stored as chars

    };
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

=head2 get_filename

	Convenience function for getting filenames of different type.

=cut

sub get_output_filename {

    (my $self, my $job_id, my $prefix, my $extension) = @_;

    return $self->{'opt'}{'temp_dir'} . '/' . $job_id . '_' . $prefix . '.' . $extension;
}

1;

