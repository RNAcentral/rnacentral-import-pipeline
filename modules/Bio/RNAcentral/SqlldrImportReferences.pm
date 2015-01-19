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

    Bio::RNAcentral::SqlldrImportReferences

=head1 DESCRIPTION

    Import literature references into the RNAcentral Oracle database.
    Since some fields (e.g. authors) are longer than 4000 characters, they are stored as clobs.
    Parallel load option is not allowed when loading lob columns (SQL*Loader-971),
    so the data are loaded sequentially.

=cut

package Bio::RNAcentral::SqlldrImportReferences;

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

=cut

sub _make_ctl_file {
    my $self = shift;

    open my $fh, '>', $self->{'local'}{'ctlfile'} or die $!;

    print $fh "LOAD DATA\n";
    $self->_list_input_files($fh, $self->{'local'}{'path'}, 'csv');

    print $fh <<CTL;
INTO TABLE $self->{'opt'}{'references_table'} TRUNCATE
FIELDS TERMINATED BY ',' enclosed by '"'
(
    MD5 char,
    ACCESSION char,
    AUTHORS char(1000000),
    LOCATION char(4000),
    TITLE char(4000),
    PMID char,
    DOI char
)
CTL
    close $fh;
}


=head2 update

    Import the data and perform pl/sql update.

=cut

sub update {
    my $self = shift;
    $self->{'logger'}->info("Loading references");
    $self->load_staging_table();
    $self->update_literature_references();
}


1;
