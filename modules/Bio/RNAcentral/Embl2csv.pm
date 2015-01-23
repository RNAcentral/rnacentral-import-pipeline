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

    This is the main program for parsing the ENA Non-coding product
    as well as other products used in the RNAcentral pipeline.

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

use Bio::SeqIO;   # BioPerl
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

    # create output directory for genomic coordinates
    mkdir $self->get_genomic_locations_path() unless -d $self->get_genomic_locations_path();

    return $self;
}


=head2 embl2csv

    Read an embl file and create the output files in csv format.
    This is the main public method of this module.

=cut

sub embl2csv {
    (my $self, my $file) = @_;

    # get file name without extension
    my $job_id = File::Basename::fileparse($file, qr/\.[^.]*/);

    my $fname_long    = $self->get_output_filename($self->get_long_folder_path(),  $job_id, 'long',  'csv');
    my $fname_short   = $self->get_output_filename($self->get_short_folder_path(), $job_id, 'short', 'csv');
    my $fname_refs    = File::Spec->catfile($self->get_refs_path(),     $job_id . '.csv');
    my $fname_ac_info = File::Spec->catfile($self->get_info_path(),     $job_id . '.csv');
    my $fname_genloc  = File::Spec->catfile($self->get_genomic_locations_path(), $job_id . '.csv');

    # open output files
    my $fh_long    = IO::File->new("> $fname_long")    or $self->{'logger'}->logdie("Couldn't open file $fname_long");
    my $fh_short   = IO::File->new("> $fname_short")   or $self->{'logger'}->logdie("Couldn't open file $fname_short");
    my $fh_refs    = IO::File->new("> $fname_refs")    or $self->{'logger'}->logdie("Couldn't open file $fname_refs");
    my $fh_ac_info = IO::File->new("> $fname_ac_info") or $self->{'logger'}->logdie("Couldn't open file $fname_ac_info");
    my $fh_gen_loc = IO::File->new("> $fname_genloc")  or $self->{'logger'}->logdie("Couldn't open file $fname_genloc");

    # open file with BioPerl
    $self->{'logger'}->info("Reading $file");
    my $stream = Bio::SeqIO->new(-file => $file, -format => 'EMBL');

    my $seq;

    # loop over seq records
    while ( eval { $seq = $stream->next_seq() } ) {

        unless ( $self->_is_valid_sequence($seq->seq) ) {
            $self->{'logger'}->info("Skipping invalid record: " . $seq->display_id);
            next;
        }

        my $xrefs = _get_xrefs($seq);

        # skip Rfam mismatches to avoid creating duplicates
        if ( skip_rfam_mismatch($seq, $xrefs) ) {
            $self->{'logger'}->info("Skipping mismatching Rfam record" . $seq->display_id);
            next;
        }

        _parse_sequences_and_xrefs($seq, $xrefs, $self->{'opt'}{'maxseqshort'}, $fh_long, $fh_short);
        _parse_literature_references($seq, $xrefs, $fh_refs);
        _parse_accession_data($seq, $xrefs, $fh_ac_info);
        _parse_genomic_locations($seq, $xrefs, $fh_gen_loc);
    }

    # if there are errors, delete already created files to prevent data import
    if ( $@ ) {
        $self->{'logger'}->logwarn('BioPerl error');
        $self->{'logger'}->logwarn($@);
        $self->{'logger'}->logwarn("Not safe to continue, skipping $file");
        unlink $fname_short, $fname_long, $fname_refs;
        return ();
    }

    $fh_short->close;
    $fh_long->close;
    $fh_refs->close;
    $fh_ac_info->close;
    $fh_gen_loc->close;

    # remove empty files
    _delete_empty_file($fname_refs);
    _delete_empty_file($fname_ac_info);

    # return an array with csv files
    return (_delete_empty_file($fname_short), _delete_empty_file($fname_long));
}


=head2 skip_rfam_mismatch

    Return 0 if a sequence should be skipped, otherwise return 1.

    Some sequences match more than one Rfam model, for example
    a bacterial SSU sequence can match up to four Rfam models:

    * RF01959 SSU_rRNA_archaea
    * RF00177 SSU_rRNA_bacteria
    * RF01960 SSU_rRNA_eukarya  Eukaryota and not Microsporidia
    * RF02542 SSU_rRNA_microsporidia Microsporidia

    When the taxon of the sequence and the Rfam model do not match,
    such entries should be skipped, so for example
    only bacterial sequences matching bacterial Rfam model are kept.

=cut

sub skip_rfam_mismatch {

    my ($seq, $xrefs) = @_;
    my $skip_entry = 0;

    if ( _is_rfam_entry($seq->display_id) ) {
        my $classification = join('; ', reverse $seq->species->classification);

        if ($$xrefs[0]->{'primary_id'}      eq 'RF01959' && $classification !~ m/Archaea/) {
            # SSU_rRNA_archaea not matching Archaea
            $skip_entry = 1;
        } elsif ($$xrefs[0]->{'primary_id'} eq 'RF00177' && $classification !~ m/Bacteria/) {
            # SSU_rRNA_bacteria not matching Bacteria
            $skip_entry = 1;
        } elsif ($$xrefs[0]->{'primary_id'} eq 'RF01960' && $classification !~ m/Eukaryota/) {
            # SSU_rRNA_eukarya not matching Eukarya
            $skip_entry = 1;
        } elsif ($$xrefs[0]->{'primary_id'} eq 'RF01960' && $classification =~ m/Microsporidia/) {
            # SSU_rRNA_eukarya matching Microsporidia
            $skip_entry = 1;
        } elsif ($$xrefs[0]->{'primary_id'} eq 'RF02542' && $classification !~ m/Microsporidia/) {
            # SSU_rRNA_microsporidia not matching Microsporidia
            $skip_entry = 1;
        }
    }

    return $skip_entry;
}


=head2 _is_valid_sequence

    Return 0 if a sequence is not valid, otherwise return 1.

    Invalid sequences:
    * shorter than minimum
    * contain only Ns
    * have high N-content

    The cutoff values are specified in Bio::RNAcentral::Base.

=cut

sub _is_valid_sequence {

    my ($self, $sequence) = @_;
    my $status;

    if (length($sequence) == 0) {
        return 0; # immediately skip empty sequences
    }

    my $N_length = () = $sequence =~ /N/ig;
    my $N_percentage = ($N_length * 100)/length($sequence);

    if (length($sequence) < $self->{'opt'}{'minseqlength'}) {
        $status = 0;
    } elsif ($N_percentage > $self->{'opt'}{'maxNcontent'} ) {
        $status = 0;
    } elsif ($N_length == length($sequence)) {
        $status = 0;
    } else {
        $status = 1;
    }

    return $status;
}


=head2 _parse_accession_data

    Prepare the data for the load_rnc_accessions table, which will be merged into
    the main rnc_accessions table.

=cut

sub _parse_accession_data {
    my ($seq, $xrefs, $fh_ac_info) = @_;
    my ($text, $accession, $parent_accession, $seq_version, $feature_location_start,
        $feature_location_end, $feature_name, $ordinal, $species,
        $common_name, $keywords, $project, $is_composite, $non_coding_id,
        $source_accession, %feature_tags, $db_xrefs);

    $source_accession = ${$xrefs}[0]{'accession'}; # accession of the source ENA entry

    $species = $seq->species;

    $project = _get_project_id($seq);

    $keywords = _sanitize($seq->keywords);
    $keywords =~ s/\.$//;

    # get feature table qualifiers
    %feature_tags = _get_feature_tags($seq);

    # get db_xref feature qualifiers
    $db_xrefs = _get_db_xref_feature_qualifiers($seq);

    # iterate over all xrefs from the entry
    for my $xref (@{$xrefs}) {

        $accession = $xref->{'accession'};

        # GU187164.1:15402..15468:tRNA
        # <accession>.<sequence version>:<feature location string>:<feature name>[:ordinal]
        $accession =~ /(\w+)\.(\d+):(\d+)\.\.(\d+):(\w+)(:(\d+))*/;
        $parent_accession       = $1;
        $seq_version            = $2;
        $feature_location_start = $3;
        $feature_location_end   = $4;
        $feature_name           = $5;
        $ordinal                = $7;

        # for RFAM entries, feature_name in the id is 'rfam', so retrieve
        # feature_name from feature table
        # same for RDP
        if ( _is_rfam_entry($accession) or _is_rdp_entry($accession) ) {
            for my $feat_object ($seq->get_SeqFeatures) {
                if ($feat_object->primary_tag eq 'source') {
                    next;
                } else {
                    $feature_name = $feat_object->primary_tag;
                }
            }
        }

        if ( _is_xref_from_dr_line($xref->{'database'}) ) {
            $is_composite = 'Y';
            $non_coding_id = $source_accession;
        } else {
            $is_composite = 'N';
            $non_coding_id = '';
        }

        $text = '"' . join('","', ($accession,
                                   $parent_accession,
                                   $seq_version,
                                   $feature_location_start,
                                   $feature_location_end,
                                   $feature_name,
                                   _sanitize($ordinal),
                                   $is_composite,
                                   $non_coding_id, # source ENA id
                                   $xref->{'database'},
                                   $xref->{'primary_id'},
                                   $xref->{'optional_id'},
                                   $project,
                                   _sanitize($seq->division),
                                   $keywords,
                                   _get_description_line($seq->desc),
                                   $species->binomial('FULL'),
                                   _sanitize($species->common_name),
                                   _sanitize($species->organelle),
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
                                   $feature_tags{'standard_name'},
                                   $db_xrefs,
                                   ) ) . "\"\n";
        print $fh_ac_info $text;
    } # end for loop
}


=head2

    Every sequence entry is an xref that is tracked in the RNAcentral system.
    In addition, the DR lines from the Expert Databases create additional xrefs
    that are also tracked.

=cut

sub _is_xref_from_dr_line {

    my $database = shift;

    my %non_dr_line_databases = (
        'ENA'    => 1,
        'RDP'    => 1,
        'RFAM'   => 1,
        'REFSEQ' => 1,
    );

    if ( exists($non_dr_line_databases{$database}) ) {
        return 0;
    } else {
        return 1;
    }
}


=head2 _get_description_line

    Tweak the description line to remove the "TPA: " string.

=cut

sub _get_description_line {

    my $description = shift;
    $description =~ s/^TPA\: //;
    return _sanitize($description);
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


=head2 _get_db_xref_feature_qualifiers

    Get the /db_xref feature qualifiers that contain feature-level xrefs
    to other resources.

    Examples of useful db_xrefs:
    * GeneID for RefSeq
    * HGNC ids for gtRNAdb

    The `source` feature is not analyzed because it only contains the taxid db_xref.

=cut

sub _get_db_xref_feature_qualifiers {

    my $seq = shift;
    my @db_xrefs;

    # loop over features
    for my $feat_object ($seq->get_SeqFeatures) {
        if ( $feat_object->primary_tag ne 'source' ) {
            # loop over feature qualifiers (tags)
            for my $tag ($feat_object->get_all_tags) {
                for my $value ($feat_object->get_tag_values($tag)) {
                    if ($tag eq 'db_xref') {
                        push @db_xrefs, $value;
                    }
                }
            }
        }
    }

    return join(' ', @db_xrefs);
}


=head2 _get_taxid

    Get taxonomic id which is found in the `db_xref` field under `source`.

=cut

sub _get_taxid {

    my $seq = shift;
    my $taxid = '';

    # loop over features
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

    return $taxid;
}

=head2 _parse_sequences_and_xrefs

    Prepare data for the load_rnacentral_all table.

=cut

sub _parse_sequences_and_xrefs {

    my ($seq, $xrefs, $maxseqshort, $fh_long, $fh_short) = @_;

    my $text     = '';
    my $sequence = $seq->seq;
    my $crc64    = SWISS::CRC64::crc64($sequence);
    my $md5      = md5_hex($sequence);

    # treat DRs as independent xrefs
    foreach my $xref ( @$xrefs ) {
        $text .=  join(',', ($crc64,
                             $seq->length,
                             $sequence,
                             $xref->{'database'},
                             $xref->{'accession'},
                             $xref->{'optional_id'},
                             $seq->seq_version,
                             _get_taxid($seq),
                             $md5)) . "\n";
    }

    # print data in different files depending on sequence length
    if ($seq->length > $maxseqshort) {
        print $fh_long $text;
    } else {
        print $fh_short $text;
    }
}


=head2 _is_rdp_entry

    Return 1 if an accession is an RDP entry, return 0 otherwise.

=cut

sub _is_rdp_entry {
    my $accession = shift;

    if ( $accession =~ /RDP$/i ) {
        return 1;
    } else {
        return 0;
    }
}

=head2 _is_rfam_entry

    Return 1 if an accession is an RFAM entry, return 0 otherwise.

=cut

sub _is_rfam_entry {
    my $accession = shift;

    if ( $accession =~ /rfam$/i ) {
        return 1;
    } else {
        return 0;
    }
}

=head _is_refseq_entry

    Return 1 if an accession is a RefSeq entry, return 0 otherwise.

=cut

sub _is_refseq_entry {
    my $keywords = shift;

    if ( $keywords =~ 'RefSeq;' ) {
        return 1;
    } else {
        return 0;
    }
}


=head2 _get_xrefs

    Read DR lines with database cross-references into an array of hashes.

=cut

sub _get_xrefs {

    my $seq = shift;

    my (@data, $database, $primary_id, $optional_id, $project, $composite_id);

    # first element is the source ENA entry unless the entry is from RFAM/RefSeq/RDP
    if ( !_is_rfam_entry($seq->display_id) and !_is_refseq_entry($seq->keywords) and !_is_rdp_entry($seq->display_id)) {
        $data[0] = { database    => 'ENA',
                     accession   => $seq->display_id,
                     primary_id  => '',
                     optional_id => '' };
    }

    # todo: remove when DR lines are added
    # set a flag if DR lines were found
    my %flags = (
        'mirbase_DR' => 0,
        'vega_DR' => 0
    );

    # append other xrefs from DR lines
    my $anno_collection = $seq->annotation;
    for my $key ( $anno_collection->get_all_annotation_keys ) {
        my @annotations = $anno_collection->get_Annotations($key);
        for my $value ( @annotations ) {
            if ( $value->tagname eq "dblink" ) {
                $database    = _sanitize($value->database());

                # skip these DR lines
                if ($database eq 'MD5' or $database =~ /SILVA/i or $database =~ /BioSample/i ) {
                    next;
                }

                # skip RFAM DR lines in regular Non-coding entries
                if (!_is_rfam_entry($seq->display_id) and $database eq 'RFAM') {
                    next;
                }

                # TODO: remove this statement when DR lines are updated
                if ($database =~ /Vega/i) {
                    # set a flag that Vega DR lines were injected
                    $flags{'vega_DR'} = 1;
                } elsif ($database =~ /mirbase/i) {
                    # set a flag that mirbase DR lines were injected
                    $flags{'mirbase_DR'} = 1;
                }

                $primary_id  = _sanitize($value->primary_id());
                $optional_id = _sanitize($value->optional_id());

                # use a shorter label for tmRNA-Website
                if ($database eq 'tmRNA-Website') {
                    $database = 'tmRNA_Web';
                }

                # RFAM/RefSeq/RDP entries already have a unique id
                if ($database eq 'RFAM' or $database eq 'RefSeq' or $database eq 'RDP') {
                    $composite_id = $seq->display_id;
                } else {
                    $composite_id = _get_expert_db_id($seq->display_id, $database, $primary_id);
                }

                push @data, {
                              accession   => $composite_id,
                              primary_id  => $primary_id, # external id
                              optional_id => $optional_id,
                              database    => uc($database), # convert to upper case
                            };
            }
        }
    }

    # todo: remove this temporary fix when DR lines are added to all entries
    # inject CC lines only if no DR lines were found
    $project = _get_project_id($seq);
    if ($flags{'mirbase_DR'} == 0 && $project eq 'PRJEB4451') {
        push @data, _inject_xrefs_manually($seq, 'MIRBASE');
    } elsif ($flags{'vega_DR'} == 0 && $project eq 'PRJEB4568') {
        push @data, _inject_vega_xrefs_manually($seq, 'Vega');
    }

    # VEGA xrefs require special treatment
    @data = @{_combine_vega_xrefs(\@data)};

    return \@data;
}


=head2 _inject_xrefs_manually

    Some xrefs are added not as DR lines but as structured comments
    because the ENA xref pipeline cannot add xrefs between ENA releases.

    As a temporary measure, these xrefs are identified using project ids
    and the DR data are parsed from the CC lines.

    Some entries from an expert database may already have DR lines while others
    may still have CC lines, depending on the status of the ENA xref pipeline.

    Example CC lines:
    lncRNAdb; 190; 7SK.
    gtRNAdb; Cand_Meth_boonei_6A8/Cand_Meth_boonei_6A8-summary.
    miRBase; MI0000182.
=cut

sub _inject_xrefs_manually {
    my ($seq, $db_name) = @_;
    my ($db_project_id, $external_id, $optional_id);

    my %project_ids = (
        'GTRNADB'  => 'PRJEB5173',
        'LNCRNADB' => 'PRJEB6238',
        'MIRBASE'  => 'PRJEB4451',
    );

    if ( exists($project_ids{$db_name} ) ) {
        $db_project_id = $project_ids{$db_name};
    }

    my $entry_project_id = _get_project_id($seq);

    if ($entry_project_id eq $db_project_id) {
        my @annotations = $seq->annotation->get_Annotations('comment');
        for my $value ( @annotations ) {
            my $comment = $value->display_text;
            my $search_string = "$db_name; (\\S+)(; (\\S+))?\\.";
            if ($comment =~ /$search_string/i) {
                $external_id = $1;
                if (defined $3) {
                    $optional_id = $3;
                } else {
                    $optional_id = '';
                }
                return {
                    accession   => _get_expert_db_id($seq->display_id, $db_name, $external_id),
                    primary_id  => $external_id,
                    optional_id => $optional_id,
                    database    => $db_name,
                };
            }
        }
    }
    return ();
}

=head2 _inject_vega_xrefs_manually

    Some Vega xrefs are found in CC lines. They are called "Vega" and not
    "Vega-Gn" or "Vega-Tr". Also, there are two lines to be extracted from
    a comment.
    TODO: revise this code once all Vega xrefs have DR lines.

    Example:
    CC   Vega; OTTHUMT00000475966; AC037193.1.
    CC   Vega; OTTHUMG00000155019; AC037193.1.

=cut

sub _inject_vega_xrefs_manually {
    my ($seq, $db_name) = @_;
    my ($db_project_id, $external_id, $optional_id, @lines);
    my @vega_xrefs = ();

    my %project_ids = (
        'Vega' => 'PRJEB4568',
    );

    if ( exists($project_ids{$db_name} ) ) {
        $db_project_id = $project_ids{$db_name};
    }

    my $entry_project_id = _get_project_id($seq);

    if ($entry_project_id eq $db_project_id) {
        my @annotations = $seq->annotation->get_Annotations('comment');
        for my $value ( @annotations ) {
            my $comment = $value->display_text;

            my @regexes = (qr/Vega; (OTTHUMG\w+);/, qr/Vega; (OTTHUMT\w+);/);

            for my $regex (@regexes) {
                if ($comment =~ $regex) {
                    $external_id = $1;
                    if (defined $3) {
                        $optional_id = $3;
                    } else {
                        $optional_id = '';
                    }
                    if ($external_id =~ /^OTTHUMG/) {
                        $db_name = 'VEGA-Gn';
                    } elsif ($external_id =~ /^OTTHUMT/) {
                        $db_name = 'VEGA-Tr';
                    }
                    push @vega_xrefs, {
                        accession   => _get_expert_db_id($seq->display_id, $db_name, $external_id),
                        primary_id  => $external_id,
                        optional_id => $optional_id,
                        database    => $db_name,
                    };
                }
            }
        }
    }

    return @vega_xrefs;
}

=head2 _get_expert_db_id

    Create unique ids for external databases.

    Format: <ENA_accession> : uppercase(<database>) : <external_id>
    The id is truncated at 100 characters because a small number of accessions
    exceeded that length, and it's not clear whether it's useful to store
    longer ids.

    Example: HG497133.1:1..472:ncRNA:VEGA:OTTHUMG00000161051
=cut

sub _get_expert_db_id {

    my ($ena_source, $database, $external_id) = @_;

    my $composite_id = $ena_source . ':' . uc($database) . ':' . $external_id;
    $composite_id =~ s/VEGA-Gn|VEGA-Tr/VEGA/i;

    my $MAX_ACCESSION_LENGTH = 100;
    return substr($composite_id, 0, $MAX_ACCESSION_LENGTH);
}


=head2 _combine_vega_xrefs

    The VEGA xrefs is defined in 2 DR lines:
    DR   VEGA-Gn; OTTHUMG00000013241; Homo sapiens. # genes
    DR   VEGA-Tr; OTTHUMT00000037008; Homo sapiens. # transcripts

    This procedure combines them as if they were in one:
    DR   VEGA; OTTHUMG00000013241; OTTHUMT00000037008.

    The species is stored elsewhere so this information is not recorded.

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


=head2 _parse_genomic_locations

    Find genomic locations relative to the primary accessions.
    These data can be used in conjunction with Ensembl Perl API to find
    toplevel genomic coordinates.

    TPA and Non-TPA entries are treated differently because TPAs refer to
    existing primary accessions and Non-TPA entries are primary accessions
    in their own right.

=cut

sub _parse_genomic_locations {

    my ($seq, $xrefs, $fh_gen_loc) = @_;
    my $text = '';
    my @locations;

    # determine if the entry is a TPA or a primary accession
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
                for my $xref (@{$xrefs}) {
                    $text .= '"' . join('","', ($xref->{'accession'},
                                                $exon->{'primary_identifier'},
                                                $exon->{'primary_start'},
                                                $exon->{'primary_end'},
                                                $exon->{'strand'})
                                        ) . "\"\n";
                }
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
                    for my $xref (@{$xrefs}) {
                        $text .= '"' . join('","', ($xref->{'accession'}, # accession used by RNAcentral
                                                    $seq->accession,      # primary identifier, e.g. AB446294.1
                                                    $location->start,     # local_start (relative to the parent accession)
                                                    $location->end,       # local_end (relative to the parent accession)
                                                    $location->strand)    # 1 or -1
                                            ) . "\"\n";
                    }
                }
            }
        }
    }
    print $fh_gen_loc $text;
}


=head2 _parse_literature_references

    Prepare data for the load_rnc_references table.

    Get csv text with literature references given a BioPerl seq object.

=cut

sub _parse_literature_references {

    my ($seq, $xrefs, $fh_refs) = @_;

    my ($authors, $location, $title, $pmid, $doi, $md5);
    my $text = '';

    my $anno_collection = $seq->annotation;
    for my $key ( $anno_collection->get_all_annotation_keys ) {
        my @annotations = $anno_collection->get_Annotations($key);
        for my $value ( @annotations ) {
            # all reference lines: RN, RC, RP, RX, RG, RA, RT, RL
            if ( $value->tagname eq 'reference' ) {
                $authors  = _sanitize($value->authors());   # RA line
                $location = _sanitize($value->location());  # RL line
                $title    = _sanitize($value->title());     # RT line
                $pmid     = _sanitize($value->pubmed());    # parsed RX line
                $doi      = _sanitize($value->doi());       # parsed RX line
                $md5      = md5_hex($authors . $location . $title);

                # iterate over all xrefs from the entry
                for my $xref (@{$xrefs}) {
                    $text .= '"' . join('","', ($md5,
                                                $xref->{'accession'},
                                                $authors,
                                                $location,
                                                $title,
                                                $pmid,
                                                $doi) ) . "\"\n";
                }
            }
        }
    }
    print $fh_refs $text;
}

####################
#                  #
# Helper functions #
#                  #
####################

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


=head2 _delete_empty_file

    Return the filename if the file is not empty, otherwise delete the file
    and return an empty list.

=cut

sub _delete_empty_file{
    my $file = shift;

    if ( -s $file > 0 ) {
        return $file;
    } else {
        unlink $file;
        return ();
    }
}


=head2 _escape_quotes

    Prepare text for writing out in comma-separated format with each field
    wrapped in quotes.

=cut

sub _escape_quotes {
    my $text = shift;

    $text =~ s/^"+//; # remove leading and trailing quotes
    $text =~ s/"+$//;

    $text =~ s/"/""/g; # escape the remaining double quotes

    return $text;
}

=head2 _sanitize

    Get a non-empty value with escaped quotes.
=cut

sub _sanitize {
    my $value = shift;
    return ( defined( $value ) ? _escape_quotes($value) : '' );
}


1;
