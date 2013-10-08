=pod

=head1 NAME

    Bio::RNAcentral::Embl2csv

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 CONTACT


=cut

package Bio::RNAcentral::Embl2csv;

use strict;
use warnings;

# override default modules, use embl.pm and CRC64.pm modules from the lib directory
use Cwd            ();
use File::Basename ();
use File::Spec;
BEGIN {
    my $root_dir = File::Basename::dirname( File::Basename::dirname( Cwd::realpath($0) ) );
    unshift @INC, "$root_dir/lib";
}

use File::Spec;
use File::Find qw(finddepth);

use Bio::SeqIO;   # BioPerl is used for reading embl files
use SWISS::CRC64; # cyclic redundancy check
use Digest::MD5 qw(md5_hex);

use base ('Bio::RNAcentral::Base');


sub new {
    my ($class, $opt) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt);

    # create output directory for references if it doesn't exist
    mkdir $self->get_refs_path() unless -d $self->get_refs_path();

    # create output directory for long files
    mkdir $self->get_long_folder_path() unless -d $self->get_long_folder_path();

    # create output directory for short files
    mkdir $self->get_short_folder_path() unless -d $self->get_short_folder_path();

    # create output directory for info files
    mkdir $self->get_info_path() unless -d $self->get_info_path();

    # create output directory for composite id files
    mkdir $self->get_comp_id_path() unless -d $self->get_comp_id_path();

    return $self;
}


=head2 embl2csv

    Read an embl file and create two csv files,
    one with short, one with long sequences.
    In addition, create csv files with database cross-references
    and literature references.

=cut

sub embl2csv {
    (my $self, my $file) = @_;

    # get file name without extension
    my $job_id = File::Basename::fileparse($file, qr/\.[^.]*/);

    my $fname_long    = $self->get_output_filename($self->get_long_folder_path(),  $job_id, 'long',  'csv');
    my $fname_short   = $self->get_output_filename($self->get_short_folder_path(), $job_id, 'short', 'csv');
    my $fname_refs    = File::Spec->catfile($self->get_refs_path(),    $job_id . '.csv');
    my $fname_ac_info = File::Spec->catfile($self->get_info_path(),    $job_id . '.csv');
    my $fname_comp_id = File::Spec->catfile($self->get_comp_id_path(), $job_id . '.csv');

    # open output files
    my $fh_long    = IO::File->new("> $fname_long")    or $self->{'logger'}->logdie("Couldn't open file $fname_long");
    my $fh_short   = IO::File->new("> $fname_short")   or $self->{'logger'}->logdie("Couldn't open file $fname_short");
    my $fh_refs    = IO::File->new("> $fname_refs")    or $self->{'logger'}->logdie("Couldn't open file $fname_refs");
    my $fh_ac_info = IO::File->new("> $fname_ac_info") or $self->{'logger'}->logdie("Couldn't open file $fname_ac_info");
    my $fh_comp_id = IO::File->new("> $fname_comp_id") or $self->{'logger'}->logdie("Couldn't open file $fname_comp_id");

    # initializations
    my ($i, $records, $seqs_long, $seqs_short, $dblinks_num, $refs_num, $data,
        $seq, $refs, $md5, $info, @dblinks);
    $i = $records = $seqs_long = $seqs_short = $dblinks_num = $refs_num = 0;
    $refs = '';

    # open file with BioPerl
    $self->{'logger'}->info("Reading $file");
    my $stream = Bio::SeqIO->new(-file => $file, -format => 'EMBL');

    # loop over seq records
    while ( eval { $seq = $stream->next_seq() } ) {
        $i++;

        $md5 = md5_hex($seq->seq);

        # get basic data for UniParc-style functionality
        $data = $self->_get_basic_data($seq, $md5);

        unless ( $data->{'isValid'} ) {
            $self->{'logger'}->logwarn("Skipping record $i: $data->{'text'}");
            next;
        }
        $dblinks_num += $data->{'dblinks_num'};

        # print the data in different files depending on sequence length
        if ( $data->{'isLong'} ) {
            print $fh_long $data->{'text'};
            $seqs_long++;
        } else {
            print $fh_short $data->{'text'};
            $seqs_short++;
        }

        # get literature references
        $refs = $self->_get_references($seq, $md5);
        if ( $refs ) {
            print $fh_refs $refs->{'text'};
            $refs_num += $refs->{'num'};
        }

        # get general info about accessions
        $info = $self->_get_ac_info($seq);
        if ( $info ) {
            print $fh_ac_info $info->{'text'};
        }

        # write out composite id data
        if ( $data->{'comp_id_text'} ) {
            print $fh_comp_id $data->{'comp_id_text'};
        }

    }

    # if there are any problems, delete already created files to prevent any import
    if ( $@ ) {
        $self->{'logger'}->logwarn('BioPerl error');
        $self->{'logger'}->logwarn($@);
        $self->{'logger'}->logwarn("Not safe to continue, skipping $file");
        unlink $fname_short, $fname_long, $fname_refs;
        return ();
    }

    $records += $i;

    $self->{'logger'}->info("$records records, $seqs_short short, $seqs_long long, " .
                            "$dblinks_num database links, $refs_num literature references");

    $fh_short->close;
    $fh_long->close;
    $fh_refs->close;
    $fh_ac_info->close;
    $fh_comp_id->close;

    $self->_check_temp_files($fname_refs);

    # return an array with csv files
    return ($self->_check_temp_files($fname_short), $self->_check_temp_files($fname_long));
}


# delete temp files if empty, return only files with content
sub _check_temp_files{
    (my $self, my $file) = @_;

    if ( -s $file > 0 ) {
        return $file;
    } else {
        $self->{'logger'}->info("File $file is empty and will be deleted");
        unlink $file;
        return ();
    }
}


=head _get_ac_info

    Get general info about accessions.

=cut

sub _get_ac_info {

    (my $self, my $seq)  = @_;

    my ($text, $accession, $seq_version, $feature_location_start,
        $feature_location_end, $feature_name, $ordinal, $species,
        $keywords, $project);

    # GU187164.1:15402..15468:tRNA
    # <accession>.<sequence version>:<feature location string>:<feature name>[:ordinal]
    $seq->display_id =~ /(\w+)\.(\d+):(\d+)\.\.(\d+):(\w+)(:(\w+))*/;
    $accession              = $1;
    $seq_version            = $2;
    $feature_location_start = $3;
    $feature_location_end   = $4;
    $feature_name           = $5;
    $ordinal                = $7;

    $species = $seq->species;

    $project = _get_project_id($seq);

    $keywords = _nvl($seq->keywords);
    $keywords =~ s/\.$//;

    $text = '"' . join('","', (_nvl($seq->display_id),
                               $accession,
                               $seq_version,
                               $feature_location_start,
                               $feature_location_end,
                               $feature_name,
                               _nvl($ordinal),
                               $project,
                               _nvl($seq->division),
                               $keywords,
                               _nvl(_sanitize($seq->desc)),
                               $species->binomial,
                               _nvl($species->organelle),
                               join('; ', reverse $species->classification),
                               ) ) . "\"\n";
    return { text => $text };
}


=head2 _get_project_id

    A wrapper to retrieve project id.

=cut

sub _get_project_id {
    my $seq = shift;

    my $project = ($seq->annotation->get_Annotations('project_id'))[0];
    if ( $project ) {
        return $project->display_text();
    } else {
        return '';
    }
}


sub _get_dblinks {

    (my $self, my $seq)  = @_;

    my (@data, $database, $primary_id, $optional_id, $project);

    # initialize as an ENA entry
    $data[0] = { database    => 'ENA',
                 primary_id  => $seq->display_id,
                 accession   => $seq->display_id,
                 optional_id => '' };

    # get mirbase entries by project ids
    # todo: remove this temporary fix when DR lines are added to mirbase entries
    $project = _get_project_id($seq);
    if ($project eq 'PRJEB4451') {
        push @data, {
            primary_id  => 'MIRBASE' . '_' . $seq->display_id, # unique xref
            accession   => 'MIRBASE',
            optional_id => '',
            database    => 'MIRBASE',
        };
    }

    # add any DR entries
    my $anno_collection = $seq->annotation;
    for my $key ( $anno_collection->get_all_annotation_keys ) {
        my @annotations = $anno_collection->get_Annotations($key);
        for my $value ( @annotations ) {
            if ( $value->tagname eq "dblink" ) {
                $database    = _nvl($value->database());
                $primary_id  = _nvl($value->primary_id());
                $optional_id = _nvl($value->optional_id());

                push @data, { # create unique ids for the DR links
                              # example: RF01271_BK006945.2:460712..467569:rRNA
                              primary_id  => $primary_id . '_' . $seq->display_id,
                              accession   => $primary_id,
                              optional_id => $optional_id,
                              database    => $database,
                            };
            }
        }
    }

    return \@data;
}


=head2 _get_basic_data

    Get minimum data required for UniParc-style functionality given a BioPerl seq object.

=cut

sub _get_basic_data {

    (my $self, my $seq, my $md5)  = @_;

    my ($ac, $length, $version, $isLong, $crc64, $taxid,
        $data, $sequence, $isValid, $text, $comp_id_text, $dblink);

    $ac = $crc64 = $taxid = $data = $sequence = $text = $comp_id_text = '';
    $length = $version = $isValid = 0;

    # get new data
    $ac       = $seq->display_id;
    $length   = $seq->length;
    $version  = $seq->seq_version;
    $sequence = $seq->seq;
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

    # quality control
    if ( $ac eq '' or $sequence eq '' or $length == 0 ) {
        $isValid = 0;
    } else {
        $isValid = 1;
    }

    if ($length > $self->{'opt'}{'maxseqshort'}) {
        $isLong = 1;
    } else {
        $isLong = 0;
    }

    # get database links
    my $dblinks = $self->_get_dblinks($seq);

    # treat DRs as independent xrefs
    foreach my $dblink ( @$dblinks ) {
        $text .=  join(',', ($crc64,
                             $length,
                             $sequence,
                             $dblink->{'database'},
                             $dblink->{'primary_id'},
                             $dblink->{'optional_id'}, # TODO: drop this column
                             $version,
                             $taxid,
                             $md5)) . "\n";
    }

    # write out composite id correspondences
    if ( scalar @$dblinks > 1 ) {
        # skip the first entry (non-composite ENA id)
        for ( my $i=1; $i < scalar @$dblinks; $i++ ) {
            $dblink = @$dblinks[$i];
            $comp_id_text .= '"' . join('","', ($dblink->{'primary_id'},  # composite id like RF01271_BK006945.2:460712..467569:rRNA
                                                $seq->display_id,         # ENA id BK006945.2:460712..467569:rRNA
                                                $dblink->{'database'}),   # e.g. RFAM
                                                $dblink->{'optional_id'}, # e.g. 5.8S_RNA
                                                $dblink->{'accession'},   # e.g. RF01271

                                        ) . "\"\n";;
        }
    } else {
        $comp_id_text = '';
    }

    return {
        'text'         => $text,
        'isLong'       => $isLong,
        'isValid'      => $isValid,
        'comp_id_text' => $comp_id_text,
        'dblinks_num'  => scalar @$dblinks,
    };
}


=head2 _get_references

    Get csv text with literature references given a BioPerl seq object.

=cut

sub _get_references {

    (my $self, my $seq, my $md5)  = @_;

    my ($ref_authors, $ref_location, $ref_title,
        $ref_pubmed,  $ref_doi,      $ref_publisher,
        $ref_editors, $ref_consortium, $md5_authors);
    my $text = '';
    my $num = 0;

    my $anno_collection = $seq->annotation;
    for my $key ( $anno_collection->get_all_annotation_keys ) {
        my @annotations = $anno_collection->get_Annotations($key);
        for my $value ( @annotations ) {
            # all reference lines: RN, RC, RP, RX, RG, RA, RT, RL
            if ( $value->tagname eq 'reference' ) {
                $num++;
                $ref_authors   = _nvl($value->authors());   # RA line
                $md5_authors   = md5_hex($ref_authors);
                $ref_location  = _nvl($value->location());  # RL line
                $ref_title     = _nvl($value->title());     # RT line
                $ref_pubmed    = _nvl($value->pubmed());    # parsed RX line
                $ref_doi       = _nvl($value->doi());       # parsed RX line
                $ref_publisher = _nvl($value->publisher()); # parsed RL line
                $ref_editors   = _nvl($value->editors());   # parsed RL line

                $text .= '"' . join('","', ($md5,
                                            $md5_authors,
                                            $ref_authors,
                                            $ref_location,
                                            $ref_title,
                                            $ref_pubmed,
                                            $ref_doi,
                                            $ref_publisher,
                                            $ref_editors) ) . "\"\n";
            }
        }
    }
    return { text => $text, num => $num };
}


sub _sanitize {
    my $text = shift;

    $text =~ s/^"+//; # remove leading and trailing quotes
    $text =~ s/"+$//;

    $text =~ s/"/""/g; # escape the remaining double quotes

    return $text;
}

=head2 _get_dblinks

    Get database cross-references from a BioPerl seq object.

=cut


# get non-empty value, named after the NVL function in Oracle
sub _nvl {
    my $value = shift;
    return ( defined( $value ) ? _sanitize($value) : '' );
}


1;