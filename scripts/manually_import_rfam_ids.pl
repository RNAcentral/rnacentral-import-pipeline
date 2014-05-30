#!/usr/bin/env perl

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

    manually_import_rfam_ids

=head1 DESCRIPTION

    For release 1.0beta RFAM ids were not imported into the database by mistake.
    This script parses the RFAM Product files and prepares files for import
    with sqlldr.

    sqlldr control file:

    LOAD DATA
    INFILE "str '\\n'"
    APPEND
    INTO TABLE TEMP_RFAM_IDS
    FIELDS TERMINATED BY ','
    (
      ACCESSION char,
      EXTERNAL_ID char,
      OPTIONAL_ID char
    )

    sqlldr command:

    sqlldr user/password@instance data=/path/to/data/file control=/path/to/control/file

=cut


use Cwd            ();
use File::Basename ();
BEGIN {
    my $root_dir = File::Basename::dirname( File::Basename::dirname( Cwd::realpath($0) ) );
    unshift @INC, "$root_dir/modules";
}

use Bio::SeqIO;


my $file = $ARGV[0];

my @fileparts = split('/', $file);

my $stream = Bio::SeqIO->new(-file => $file, -format => 'EMBL');

my $fh = IO::File->new("> " . $fileparts[-1]);

# loop over seq records
while ( eval { $seq = $stream->next_seq() } ) {
    my $anno_collection = $seq->annotation;
    for my $key ( $anno_collection->get_all_annotation_keys ) {
        my @annotations = $anno_collection->get_Annotations($key);
        for my $value ( @annotations ) {
            if ( $value->tagname eq "dblink" and $value->database() eq 'RFAM') {
                $accession = $seq->display_id;
                $primary_id  = $value->primary_id();
                $optional_id = $value->optional_id();
                print $fh "$accession,$primary_id,$optional_id\n";
            }
        }
    }
}

$fh->close;
