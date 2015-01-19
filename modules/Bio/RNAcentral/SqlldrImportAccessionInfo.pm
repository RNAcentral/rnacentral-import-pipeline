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

    Bio::RNAcentral::SqlldrImportAccessionInfo

=head1 DESCRIPTION


=cut

package Bio::RNAcentral::SqlldrImportAccessionInfo;

use strict;
use warnings;

use File::Spec;

use base ('Bio::RNAcentral::SqlldrImportBase');


sub new {
    my ($class, $opt, $path, $prefix) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt, $path, $prefix);

    return $self;
}


=head2 _make_ctl_file

    Create a control file used by sqlldr.

    Default for `char` with unspecified length is 255 characters.
=cut

sub _make_ctl_file {
    my $self = shift;

    open my $fh, '>', $self->{'local'}{'ctlfile'} or die $!;

    print $fh "LOAD DATA\n";
    $self->_list_input_files($fh, $self->{'local'}{'path'}, 'csv');

    print $fh <<CTL;
INTO TABLE $self->{'opt'}{'ac_info_table'} TRUNCATE
FIELDS TERMINATED BY ',' enclosed by '"'
(
    ACCESSION char,
    PARENT_AC char,
    SEQ_VERSION integer external,
    FEATURE_START integer external,
    FEATURE_END integer external,
    FEATURE_NAME char,
    ORDINAL integer external,
    IS_COMPOSITE,
    NON_CODING_ID,
    DATABASE,
    EXTERNAL_ID,
    OPTIONAL_ID,
    PROJECT char,
    DIVISION char,
    KEYWORDS char,
    DESCRIPTION char,
    SPECIES char,
    COMMON_NAME char,
    ORGANELLE char,
    CLASSIFICATION char(500),
    ALLELE char,
    ANTICODON char,
    CHROMOSOME char,
    EXPERIMENT char(500),
    FUNCTION char(2000),
    GENE char,
    GENE_SYNONYM char(400),
    INFERENCE char,
    LOCUS_TAG char,
    MAP char,
    MOL_TYPE char,
    NCRNA_CLASS char,
    NOTE char(1500),
    OLD_LOCUS_TAG char,
    OPERON char,
    PRODUCT char,
    PSEUDOGENE char,
    STANDARD_NAME char,
    DB_XREF char
)
CTL
    close $fh;
}


=head2 update

    Import the data and perform pl/sql update.

=cut

sub update {
    my $self = shift;
    $self->{'logger'}->info("Loading accession info");
    $self->load_staging_table();
    $self->update_accession_info();
}


1;
