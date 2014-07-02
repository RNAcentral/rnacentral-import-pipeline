=pod

=head1 LICENSE

    Copyright [2009-2014] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

=head1 NAME

    Bio::RNAcentral::Embl2csv

=head1 DESCRIPTION


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
use JSON;
use LWP::Simple;

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

    # create output directory for assembly data files
    mkdir $self->get_assembly_path() unless -d $self->get_assembly_path();

    # create output directory for genomic locations
    mkdir $self->get_genomic_locations_path() unless -d $self->get_genomic_locations_path();

    $self->{'assemblies'} = {};

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
    my $fname_refs    = File::Spec->catfile($self->get_refs_path(),     $job_id . '.csv');
    my $fname_ac_info = File::Spec->catfile($self->get_info_path(),     $job_id . '.csv');
    my $fname_comp_id = File::Spec->catfile($self->get_comp_id_path(),  $job_id . '.csv');
    my $fname_as_info = File::Spec->catfile($self->get_assembly_path(), $job_id . '.csv');
    my $fname_genloc  = File::Spec->catfile($self->get_genomic_locations_path(), $job_id . '.csv');

    # open output files
    my $fh_long    = IO::File->new("> $fname_long")    or $self->{'logger'}->logdie("Couldn't open file $fname_long");
    my $fh_short   = IO::File->new("> $fname_short")   or $self->{'logger'}->logdie("Couldn't open file $fname_short");
    my $fh_refs    = IO::File->new("> $fname_refs")    or $self->{'logger'}->logdie("Couldn't open file $fname_refs");
    my $fh_ac_info = IO::File->new("> $fname_ac_info") or $self->{'logger'}->logdie("Couldn't open file $fname_ac_info");
    my $fh_comp_id = IO::File->new("> $fname_comp_id") or $self->{'logger'}->logdie("Couldn't open file $fname_comp_id");
    my $fh_as_info = IO::File->new("> $fname_as_info") or $self->{'logger'}->logdie("Couldn't open file $fname_as_info");
    my $fh_gen_loc = IO::File->new("> $fname_genloc")  or $self->{'logger'}->logdie("Couldn't open file $fname_genloc");

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
            $self->{'logger'}->info("Skipping invalid record $i: $data->{'text'}");
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
        $refs = $self->_get_references($seq);
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

        # get genomic locations
        print $fh_gen_loc $self->_get_genomic_locations($seq);

        # assembly information with genome locations
        my $assembly_json = ($seq->annotation->get_Annotations('assembly_json'))[0];
        if ($assembly_json) {
            my $assembly_info = decode_json($assembly_json->display_text);
            for my $exon (@$assembly_info) {
                print $fh_as_info '"' . join( '","', ($seq->display_id,
                                                      $exon->{'local_start'},
                                                      $exon->{'local_end'},
                                                      $exon->{'primary_identifier'},
                                                      $exon->{'primary_start'},
                                                      $exon->{'primary_end'},
                                                      $exon->{'strand'})) . "\"\n";
            }
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
    $fh_as_info->close;
    $fh_gen_loc->close;

    # remove empty files
    $self->_check_temp_files($fname_refs);
    $self->_check_temp_files($fname_comp_id);
    $self->_check_temp_files($fname_as_info);
    $self->_check_temp_files($fname_ac_info);

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
        $common_name, $keywords, $project);

    # GU187164.1:15402..15468:tRNA
    # <accession>.<sequence version>:<feature location string>:<feature name>[:ordinal]
    $seq->display_id =~ /(\w+)\.(\d+):(\d+)\.\.(\d+):(\w+)(:(\w+))*/;
    $accession              = $1;
    $seq_version            = $2;
    $feature_location_start = $3;
    $feature_location_end   = $4;
    $feature_name           = $5;
    $ordinal                = $7;

    # for RFAM entries, feature_name in the id is 'rfam', so retrieve
    # feature_name from feature table
    if ($seq->display_id =~ /rfam$/) {
        for my $feat_object ($seq->get_SeqFeatures) {
            if ($feat_object->primary_tag eq 'source') {
                next;
            } else {
                $feature_name = $feat_object->primary_tag;
            }
        }
    }

    $species = $seq->species;

    $project = _get_project_id($seq);

    $keywords = _nvl($seq->keywords);
    $keywords =~ s/\.$//;

    # get feature table qualifiers
    my %feature_tags = _get_feature_tags($seq);

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
                               $species->binomial(),
                               _nvl($species->common_name),
                               _nvl($species->organelle),
                               join('; ', reverse $species->classification),
                               $feature_tags{'allele'},
                               $feature_tags{'anticodon'},
                               $feature_tags{'chromosome'},
                               $feature_tags{'experiment'},
                               $feature_tags{'function'},
                               $feature_tags{'gene'},
                               $feature_tags{'gene_synonym'},
                               $feature_tags{'inference'},
                               $feature_tags{'locus_tag'},
                               $feature_tags{'map'},
                               $feature_tags{'mol_type'},
                               $feature_tags{'ncRNA_class'},
                               $feature_tags{'note'},
                               $feature_tags{'old_locus_tag'},
                               $feature_tags{'operon'},
                               $feature_tags{'product'},
                               $feature_tags{'pseudogene'},
                               $feature_tags{'standard_name'}
                               ) ) . "\"\n";
    return { text => $text };
}


=head2 _get_feature_tags

    Retrieve a dictionary with feature data.

    The Non-coding product includes the following features:

    source
    tRNA
    tmRNA
    rRNA
    precursor_RNA
    ncRNA
    misc_RNA

    The list of feature qualifiers for these features is available at:
    http://www.insdc.org/documents/feature_table.html

=cut

sub _get_feature_tags {

    my $seq = shift;
    my @tags = (
        'allele',        # name of the allele for the given gene
        'anticodon',     # location of the anticodon of tRNA and the amino acid for which it codes
        'chromosome',    # chromosome (e.g. Chromosome number) from which the sequence was obtained
        'experiment',    # a brief description of the nature of the experimental evidence that supports the feature identification or assignment.
        'function',      # function attributed to a sequence
        'gene',          # symbol of the gene corresponding to a sequence region
        'gene_synonym',  # synonymous, replaced, obsolete or former gene symbol
        'inference',     # a structured description of non-experimental evidence that supports the feature identification or assignment.
        'locus_tag',     # a submitter-supplied, systematic, stable identifier for a gene and its associated features, used for tracking purposes
        'map',           # genomic map position of feature
        'mol_type',      # in vivo molecule type of sequence
        'ncRNA_class',   # INSDC http://www.insdc.org/rna_vocab.html
        'note',          # any comment or additional information
        'old_locus_tag', # feature tag assigned for tracking purposes
        'operon',        # name of the group of contiguous genes transcribed into a single transcript to which that feature belongs.
        'product',       # name of the product associated with the feature
        'pseudogene',    # indicates that this feature is a pseudogene of the element named by the feature key
        'standard_name', # accepted standard name for this feature
    );
    my %tag_values = ();

    # initialize with empty strings
    for my $tag (@tags) {
        $tag_values{$tag} = '';
    }

    # get data
    for my $tag (@tags) {
        for my $feat_object ($seq->get_SeqFeatures) {
            if ($feat_object->has_tag($tag)) {
                $tag_values{$tag} = join(' ', $feat_object->get_tag_values($tag));
            }
        }
    }

    return %tag_values;
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

=head2 _inject_xrefs

    Some entries don't have DR lines for all xrefs.
    As a temporary measure, these xrefs are identified using project ids
    and the DR data are parsed from the CC lines.

=cut

sub _inject_xrefs {
    my ($seq, $db_name) = @_;
    my ($db_project_id, $external_id, $optional_id);

    if ($db_name eq 'GTRNADB') {
        $db_project_id = 'PRJEB5173';
    } elsif ($db_name eq 'LNCRNADB') {
        $db_project_id = 'PRJEB6238';
    }

    my $entry_project_id = _get_project_id($seq);
    if ($entry_project_id eq $db_project_id) {
        #parse comment lines to extract xref information
        #lncRNAdb; 190; 7SK. Text continues.
        #
        # or
        #
        #Specialist DB   : gtRNAdb (Genomic tRNA Database)
        #URL             : http://lowelab.ucsc.edu/GtRNAdb/
        #gtRNAdb; Cand_Meth_boonei_6A8/Cand_Meth_boonei_6A8-summary.
        my @annotations = $seq->annotation->get_Annotations('comment');
        for my $value ( @annotations ) {
            my $comment = $value->display_text;
            my $search_string = "$db_name; (\\S+)(; (\\S+))?\\.";
            if ($comment =~ /$search_string/i) {;
                $external_id = $1;
                if (defined $3) {
                    $optional_id = $3;
                } else {
                    $optional_id = '';
                }
            } else {
                $external_id = $db_name;
                $optional_id = '';
            }
        }

        return {
            accession   => _get_composite_id($seq->display_id, $db_name, $external_id),
            primary_id  => $external_id,
            optional_id => $optional_id,
            database    => $db_name,
        };
    } else {
        return ();
    }
}


=head2 _get_xrefs

    Read DR lines with database cross-references from a BioPerl seq object.

=cut

sub _get_xrefs {

    (my $self, my $seq)  = @_;

    my (@data, $database, $primary_id, $optional_id, $project, $is_rfam_entry,
        $composite_id);

    if ( $seq->display_id =~ /rfam$/i ) {
        $is_rfam_entry = 1;
    } else {
        $is_rfam_entry = 0;
    }

    # first element is the source ENA entry unless it's an RFAM entry
    if ( !$is_rfam_entry ) {
        $data[0] = { database    => 'ENA',
                     primary_id  => $seq->display_id,
                     accession   => $seq->display_id,
                     optional_id => '' };
    }

    # get gtRNAdb and lncRNAdb entries by project ids
    # todo: remove this temporary fix when DR lines are added to all entries
    push @data, _inject_xrefs($seq, 'GTRNADB');
    push @data, _inject_xrefs($seq, 'LNCRNADB');

    # append other xrefs from DR lines
    my $anno_collection = $seq->annotation;
    for my $key ( $anno_collection->get_all_annotation_keys ) {
        my @annotations = $anno_collection->get_Annotations($key);
        for my $value ( @annotations ) {
            if ( $value->tagname eq "dblink" ) {
                $database    = _nvl($value->database());

                # skip these DR lines
                if ($database eq 'MD5' or $database =~ /SILVA/) {
                    next;
                }

                # skip RFAM DR lines in regular Non-coding entries
                if (!$is_rfam_entry and $database eq 'RFAM') {
                    next;
                }

                $primary_id  = _nvl($value->primary_id());
                $optional_id = _nvl($value->optional_id());

                # use a shorter label for tmRNA-Website
                if ($database eq 'tmRNA-Website') {
                    $database = 'tmRNA_Web';
                }

                # RFAM entries already have a unique id
                if ($database eq 'RFAM') {
                    $composite_id = $seq->display_id;
                } else {
                    $composite_id = _get_composite_id($seq->display_id, $database, $primary_id);
                }

                push @data, {
                              accession   => $composite_id,
                              primary_id  => $primary_id, # external id
                              optional_id => $optional_id,
                              database    => uc($database),
                            };
            }
        }
    }

    # VEGA xrefs require special treatment
    @data = @{_combine_vega_xrefs(\@data)};

    return \@data;
}


=head2 _get_composite_id

    Create unique ids for external databases.

    Format: <ENA_accession> : uppercase(<database>) : <external_id>
    The id is truncated at 100 characters because a small number of accessions
    exceeded that length, and it's not clear whether it's useful to store
    longer ids.

    Example: HG497133.1:1..472:ncRNA:VEGA:OTTHUMG00000161051
=cut

sub _get_composite_id {

    my ($ena_source, $database, $external_id) = @_;

    my $composite_id = $ena_source . ':' . uc($database) . ':' . $external_id;
    $composite_id =~ s/VEGA-Gn|VEGA-Tr/VEGA/i;

    my $MAX_ACCESSION_LENGTH = 100;
    return substr($composite_id, 0, $MAX_ACCESSION_LENGTH);
}


=head2 _combine_vega_xrefs

    VEGA xref is split into 2 DR lines:
    DR   VEGA-Gn; OTTHUMG00000013241; Homo sapiens. # genes
    DR   VEGA-Tr; OTTHUMT00000037008; Homo sapiens. # transcripts

    This procedure collates them as if they were in one:
    DR   VEGA; OTTHUMG00000013241; OTTHUMT00000037008.

    The "Homo sapiens" bit is not important because the taxon
    also appears elsewhere in the entry.

=cut

sub _combine_vega_xrefs {

    my ($data) = @_;
    my @new_data;
    my $new_vega_dr_link = {
        primary_id  => '',
        accession   => '',
        optional_id => '',
        database    => 'VEGA',
    };

    for my $dr_link (@$data) {
        if ( $dr_link->{'database'} =~ /^VEGA-Gn$/i ) {
            $new_vega_dr_link->{'primary_id'} = $dr_link->{'primary_id'};
            $new_vega_dr_link->{'accession'} = $dr_link->{'accession'};
        } elsif ( $dr_link->{'database'} =~ /^VEGA-Tr$/i ) {
            $new_vega_dr_link->{'optional_id'} = $dr_link->{'primary_id'};
        } else {
            push @new_data, $dr_link;
        }
    }

    if ( $new_vega_dr_link->{'primary_id'} ne '' ) {
        push @new_data, $new_vega_dr_link;
    }

    return \@new_data;
}


=head2 _is_valid_sequence

    Discard sequences shorter than minimum or those composed only of Ns.

=cut

sub _is_valid_sequence {

    my ($self, $sequence) = @_;
    my $status;

    my $N_characters = () = $sequence =~ /N/ig;
    my $N_content = ($N_characters * 100)/length($sequence);

    if (length($sequence) < $self->{'opt'}{'minseqlength'}) {
        $status = 0;
    } elsif ($N_content > $self->{'opt'}{'maxNcontent'} ) {
        $status = 0;
    } elsif ($N_characters == length($sequence)) {
        $status = 0;
    } else {
        $status = 1;
    }

    return $status;
}


=head2 _get_genomic_locations

    Find genomic locations relative to the primary accessions.
    These data can be used in conjunction with Ensembl Perl API to find
    toplevel genomic coordinates.

    TPA and Non-TPA entries are treated differently because TPAs refer to
    existing primary accessions and Non-TPA entries are primary accessions
    in their own right.

=cut

sub _get_genomic_locations {

    my ($self, $seq) = @_;
    my $text = '';
    my @locations;

    my $is_TPA;
    if ($seq->keywords =~ /TPA;/) {
        $is_TPA = 1;
    } else {
        $is_TPA = 0;
    }

    if ($is_TPA) {
        # look in the AS lines
        my $assembly_json = ($seq->annotation->get_Annotations('assembly_json'))[0];
        if ($assembly_json) {
            my $assembly_info = decode_json($assembly_json->display_text);
            for my $exon (@$assembly_info) {
                $text .= '"' . join('","', ($seq->display_id,
                                            $exon->{'primary_identifier'},
                                            $exon->{'primary_start'},
                                            $exon->{'primary_end'},
                                            $exon->{'strand'})
                                    ) . "\"\n";
            }
        }
    } else {
        # look in the feature table
        for my $feat_object ($seq->get_SeqFeatures) {
            # find tRNA, ncRNA and other non-coding features
            if ( $feat_object->primary_tag ne 'source' ) {
                # split location or a simple range
                if ( $feat_object->location->isa('Bio::Location::SplitLocationI') ) {
                    @locations = $feat_object->location->sub_Location;
                } else {
                    @locations = ($feat_object->location);
                }
                for my $location ( @locations ) {
                    $text .= '"' . join('","', ($seq->display_id,
                                                $seq->accession,   # primary identifier, e.g. AB446294.1
                                                $location->start,  # local_start (relative to the parent accession)
                                                $location->end,    # local_end (relative to the parent accession)
                                                $location->strand) # 1 or -1
                                        ) . "\"\n";
                }
            }
        }
    }
    return $text;
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

    # loop over features
    for my $feat_object ($seq->get_SeqFeatures) {
        # get taxid
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
    $isValid = $self->_is_valid_sequence($sequence);

    if ($length > $self->{'opt'}{'maxseqshort'}) {
        $isLong = 1;
    } else {
        $isLong = 0;
    }

    # get database links
    my $dblinks = $self->_get_xrefs($seq);

    # treat DRs as independent xrefs
    foreach my $dblink ( @$dblinks ) {
        $text .=  join(',', ($crc64,
                             $length,
                             $sequence,
                             $dblink->{'database'},
                             $dblink->{'accession'},
                             $dblink->{'optional_id'},
                             $version,
                             $taxid,
                             $md5)) . "\n";
    }

    # write out composite id correspondences
    if ( scalar @$dblinks > 1 ) {
        # skip the first entry (non-composite ENA id or RFAM product)
        for ( my $i=1; $i < scalar @$dblinks; $i++ ) {
            $dblink = @$dblinks[$i];
            $comp_id_text .= '"' . join('","', ($dblink->{'accession'},   # composite id like HG497133.1:1..472:ncRNA:VEGA:OTTHUMG00000161051
                                                $seq->display_id,         # ENA source accession BK006945.2:460712..467569:rRNA
                                                $dblink->{'database'}),   # e.g. RFAM
                                                $dblink->{'optional_id'}, # e.g. 5.8S_RNA
                                                $dblink->{'primary_id'},  # e.g. RF01271

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

    (my $self, my $seq)  = @_;

    my ($ref_authors, $ref_location, $ref_title,
        $ref_pmid,    $ref_doi,      $md5);
    my $text = '';
    my $num = 0;

    my $anno_collection = $seq->annotation;
    for my $key ( $anno_collection->get_all_annotation_keys ) {
        my @annotations = $anno_collection->get_Annotations($key);
        for my $value ( @annotations ) {
            # all reference lines: RN, RC, RP, RX, RG, RA, RT, RL
            if ( $value->tagname eq 'reference' ) {
                $num++;
                $ref_authors  = _nvl($value->authors());   # RA line
                $ref_location = _nvl($value->location());  # RL line
                $ref_title    = _nvl($value->title());     # RT line
                $ref_pmid     = _nvl($value->pubmed());    # parsed RX line
                $ref_doi      = _nvl($value->doi());       # parsed RX line
                $md5          = md5_hex($ref_authors . $ref_location . $ref_title);

                $text .= '"' . join('","', ($md5,
                                            $seq->display_id,
                                            $ref_authors,
                                            $ref_location,
                                            $ref_title,
                                            $ref_pmid,
                                            $ref_doi) ) . "\"\n";
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

# get non-empty value, named after the NVL function in Oracle
sub _nvl {
    my $value = shift;
    return ( defined( $value ) ? _sanitize_text($value) : '' );
}


1;
