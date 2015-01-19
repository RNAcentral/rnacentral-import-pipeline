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

    Bio::RNAcentral::SqlldrImportCoordinates

=head1 DESCRIPTION


=cut

package Bio::RNAcentral::SqlldrImportCoordinates;

use strict;
use warnings;

use base ('Bio::RNAcentral::SqlldrImportBase');


sub new {
    my ($class, $opt) = @_;

    my $folder = 'genome_coordinates';

    # run parent constructor
    my $self = $class->SUPER::new($opt, $folder);

    return $self;
}


=head2 load_staging_table

    Load the data into the staging table.

=cut

sub load_staging_table {
    my $self = shift;

    $self->{'logger'}->info("Loading genomic coordinates");

    # create one sqlldr control file
    $self->_make_ctl_file();

    # prepare sqlldr command
    my $cmd = $self->_get_sqlldr_command_list_files();

    # run sqlldr
    my $problems = $self->_run_sqlldr($cmd);

    # clean up if no errors in sqlldr
    unless ( $self->_errors_found() or $problems ) {
        unlink $self->{'local'}{'logfile'}, $self->{'local'}{'badfile'};
    }
}


=head2 _make_ctl_file

    Create a control file for sqlldr.

    List all csv files to be imported in the Control file.

    Default for `char` with unspecified length is 255 characters.
=cut

sub _make_ctl_file {
    my $self = shift;

    open my $fh, '>', $self->{'local'}{'ctlfile'} or die $!;

    print $fh "LOAD DATA\n";

    my $path = $self->get_genomic_locations_path();
    $self->_list_input_files($fh, $path, 'csv');

print $fh <<CTL;
INTO TABLE $self->{'opt'}{'coordinates_table'} TRUNCATE
FIELDS TERMINATED BY ',' enclosed by '"'
(
    ACCESSION char,
    PRIMARY_ACCESSION char,
    LOCAL_START integer external,
    LOCAL_END integer external,
    STRAND integer external
)
CTL
    close $fh;
}


=head2 update_ac_info

    Merge the accession information from the staging table.

=cut

sub update_coordinates {
    my $self = shift;

    $self->{'logger'}->info("Launching PLSQL coordinates update");

    my $command = <<PLSQL;
BEGIN

    -- remove old entries so that they can be replaced with new data
    EXECUTE IMMEDIATE '
    DELETE FROM rnc_coordinates t1
    WHERE EXISTS (
      SELECT t2.accession
      FROM load_rnc_coordinates t2
      WHERE t1.accession = t2.accession
    )
    ';

    -- move all data from the staging table
    EXECUTE IMMEDIATE '
    INSERT INTO rnc_coordinates
    (
        accession,
        primary_accession,
        local_start,
        local_end,
        name,
        primary_start,
        primary_end,
        strand
    )
    SELECT
        accession,
        primary_accession,
        local_start,
        local_end,
        name,
        primary_start,
        primary_end,
        strand
    FROM load_rnc_coordinates
    ';

    EXECUTE IMMEDIATE 'COMMIT';

    DBMS_OUTPUT.put_line('PLSQL: Coordinates updated');

END;
PLSQL

    $self->db_oracle_connect();
    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("Coordinates update failed: " . $DBI::errstr);
    $self->log_plsql_output();
    $self->{'logger'}->info("Coordinates update complete");
    $self->db_oracle_disconnect();
}


=head2 update

    Import the data into the staging table
    and move it into the production table.

=cut

sub update {
    my $self = shift;

    $self->load_staging_table();
    $self->update_coordinates();
}


1;
