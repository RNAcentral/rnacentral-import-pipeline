=pod

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION

=cut

package Bio::RNAcentral::Base;

use strict;

use Log::Log4perl;
use File::Spec;


=head2 new

	Constructor, initializes logging and sets default options.

=cut

sub new {
    my ($class, @args) = @_;

    my $self = bless {}, $class;

    $self->{'logger'} = initialize_logger();
    $self->{'opt'}    = default_options();

    return $self;
}


sub default_options {

    return {

        'file_extension'   => 'ncr',            # look for .ncr files
        'file_size_cutoff' => 50 * 10**6,       # file size cutoff in bytes
        'staging_table' => 'load_rnacentral2',  # staging table
        'release_table' => 'rnc_release',       # keeps track of all RNAcentral releases
		'maxseqlong'    => 1000000,  # maximum length for long sequences stored as clobs
		'maxseqshort'   => 4000,     # maximum length for short sequences stored as chars

    };
}


=head2 initialize_logger

    Initialize the logger object and return it so it can be reused.
    TODO: switch to logging to file or both to screen and file

=cut

sub initialize_logger {

    my $log_conf = q(
       log4perl.rootLogger              = DEBUG, LOG1
       log4perl.appender.LOG1           = Log::Log4perl::Appender::Screen
       log4perl.appender.Screen.stderr  = 0
       log4perl.appender.LOG1.layout    = Log::Log4perl::Layout::PatternLayout
       log4perl.appender.LOG1.layout.ConversionPattern = %d %p %m %n
    );

    Log::Log4perl::init_once(\$log_conf);

    my $logger = Log::Log4perl->get_logger();

    return $logger;
}

=head2 get_filename

	Convenience function for getting filenames of different type.

=cut

sub get_output_filename {

    my ($self, $path, $job_id, $prefix, $extension) = @_;

    return File::Spec->catfile($path, $job_id . '_' . $prefix . '.' . $extension);
}

1;

