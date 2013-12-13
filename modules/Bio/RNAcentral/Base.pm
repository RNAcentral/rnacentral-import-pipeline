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
    my ($class, $opt) = @_;

    my $self = bless {}, $class;

    $self->{'user'}          = $opt->{'user'};
    $self->{'password'}      = $opt->{'password'};
    $self->{'sid'}           = $opt->{'sid'};
    $self->{'port'}          = $opt->{'port'};
    $self->{'host'}          = $opt->{'host'};
    $self->{'output_folder'} = $opt->{'output_folder'};
    $self->{'release_type'}  = $opt->{'release_type'};

    $self->{'logger'} = initialize_logger();
    $self->{'opt'}    = default_options();

    return $self;
}


sub default_options {
    return {
        'file_extension'   => 'ncr',                     # look for .ncr files
        'file_size_cutoff' => 50 * 10**6,                # file size cutoff in bytes
        'staging_table'    => 'load_rnacentral_all',     # staging table
        'release_table'    => 'rnc_release',             # keeps track of all RNAcentral releases
        'references_table' => 'load_rnc_references',     # table with literature references
        'ac_info_table'    => 'rnc_ac_info_all',         # data about accessions
        'comp_id_table'    => 'rnc_composite_ids_all',   # composite ids derived from DR lines
    		'maxseqlong'       => 1000000,                   # maximum length for long sequences stored as clobs
    		'maxseqshort'      => 4000,                      # maximum length for short sequences stored as chars
        'ebi_ftp_site'     => 'ftp.ebi.ac.uk',           # Non-coding Product FTP site details
        'ebi_ftp_user'     => 'anonymous',
        'ebi_ftp_password' => '',
        'ebi_ftp_non_coding_product_folder'  => '/pub/databases/ena/non-coding',
        'ebi_ftp_non_coding_product_update'  => '/pub/databases/ena/non-coding/update',
        'ebi_ftp_non_coding_product_release' => '/pub/databases/ena/non-coding/release',
    };
}


=head2 initialize_logger

    Initialize the logger object and return it so it can be reused.
    TODO: switch to logging to file or both to screen and file

=cut

sub initialize_logger {

    my $log_conf = q(
       log4perl.rootLogger              = DEBUG, LOG1
       log4perl.appender.LOG1           = Log::Log4perl::Appender::File
       log4perl.appender.LOG1.filename  = log/rnacentral_import.log
       log4perl.appender.LOG1.autoflush = 1
       log4perl.appender.Screen.stderr  = 0
       log4perl.appender.LOG1.layout    = Log::Log4perl::Layout::PatternLayout
       log4perl.appender.LOG1.layout.ConversionPattern = %d %p %m %n

       log4perl.appender.Syncer            = Log::Log4perl::Appender::Synchronized
       log4perl.appender.Syncer.appender   = LOG1
    );

    Log::Log4perl::init_once(\$log_conf);

    my $logger = Log::Log4perl->get_logger();

    return $logger;
}


sub get_assembly_path {
  my $self = shift;
  return File::Spec->catfile($self->{'output_folder'}, 'assemblies');
}

sub get_sqlldr_bad_path {
  my $self = shift;
  return File::Spec->catfile($self->{'output_folder'}, 'bad');
}


sub get_sqlldr_log_path {
  my $self = shift;
  return File::Spec->catfile($self->{'output_folder'}, 'log');
}


sub get_short_folder_path {
  my $self = shift;
  return File::Spec->catfile($self->{'output_folder'}, 'short');
}


sub get_long_folder_path {
  my $self = shift;
  return File::Spec->catfile($self->{'output_folder'}, 'long');
}


sub get_chunks_path {
    my $self = shift;
    return File::Spec->catfile($self->{'output_folder'}, 'chunks');
}


=head2 get_comp_id_path

  Get full path to the folder containing files downloaded over ftp.

=cut

sub get_comp_id_path {
    my $self = shift;
    return File::Spec->catfile($self->{'output_folder'}, 'composite_ids');
}


=head2 get_ftp_downloads_path

  Get full path to the folder containing files downloaded over ftp.

=cut

sub get_ftp_downloads_path {
    my $self = shift;
    return File::Spec->catfile($self->{'output_folder'}, 'ncproduct');
}


=head2 get_refs_path

  Get full path of the folder with csv files containing literature
  references.

=cut

sub get_refs_path {
    my $self = shift;
    return File::Spec->catfile($self->{'output_folder'}, 'refs');
}


=head2 get_info_path

  Get full path of the folder with csv files containing general information
  about accessions, such as descriptions etc.

=cut

sub get_info_path {
    my $self = shift;
    return File::Spec->catfile($self->{'output_folder'}, 'ac_info');
}


=head2 get_filename

	Convenience function for getting filenames of different type.

=cut

sub get_output_filename {

    my ($self, $path, $job_id, $prefix, $extension) = @_;

    return File::Spec->catfile($path, $job_id . '_' . $prefix . '.' . $extension);
}

1;

