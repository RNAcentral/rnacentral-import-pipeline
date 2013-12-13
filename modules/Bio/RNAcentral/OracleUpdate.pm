#!/usr/bin/env perl

package Bio::RNAcentral::OracleUpdate;

=pod

=head1 NAME

    Bio::RNAcentral::OracleUpdate

=head1 DESCRIPTION

    Manage connection to the Oracle database
    Truncate staging table in the beginning of the update
    Create new release
    Launch PL/SQL update

=cut

use strict;
use warnings;

use DBI;
use DBD::Oracle;

use base ('Bio::RNAcentral::Base');


sub new {
    my ($class, $opt) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt);

    return $self;
}


=head2 db_oracle_connect

    Establish a database connection and keep it open.

=cut

sub db_oracle_connect {
    my $self = shift;

    my $dsn = "dbi:Oracle:host=$self->{'host'};service_name=$self->{'sid'};port=$self->{'port'}";
    my $dbh = DBI->connect($dsn, $self->{'user'}, $self->{'password'})
              or $self->{'logger'}->logdie( $DBI::errstr . "\n" );

    $self->{'dbh'} = $dbh;
    $self->{'dbh'}->func( 1000000, 'dbms_output_enable' );
    $self->{'logger'}->info("Connected to the database");
}

=head2 db_oracle_disconnect

    Close database connection.

=cut

sub db_oracle_disconnect {
    my $self = shift;

    if ( exists $self->{'dbh'} ) {
        $self->{'dbh'}->disconnect();
        $self->{'logger'}->info("Disconnected from the database");
    }
}


=head2 truncate_table

    Generic truncate table function.

=cut

sub truncate_table {
    (my $self, my $table) = @_;

    my $command = "TRUNCATE TABLE $table";
    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("Couldn't truncate $table");

    $self->{'logger'}->info("Table $table truncated");
}


=head2 update_ac_info

    Merge the accession information from the staging table.

=cut

sub update_accession_info {
    my $self = shift;

    $self->{'logger'}->info("Launching accession info update");

    my $command = <<PLSQL;
BEGIN
  RNC_UPDATE.update_accession_info();
END;
PLSQL

    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("Accession info update failed " . $DBI::errstr);

    $self->log_plsql_output();
    $self->{'logger'}->info("Accession info update complete");
}


=head2 update_composite_ids

    Merge the composite id data from the staging table.

=cut

sub update_composite_ids {
    my $self = shift;

    $self->{'logger'}->info("Launching composite id update");

    my $command = <<PLSQL;
BEGIN
  RNC_UPDATE.update_composite_ids();
END;
PLSQL

    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("Composite id update failed " . $DBI::errstr);

    $self->log_plsql_output();
    $self->{'logger'}->info("Composite id update complete");
}


=head2 update_literature_references

    Merge the literature references data from the staging table.

=cut

sub update_literature_references {
    my $self = shift;

    $self->{'logger'}->info("Launching literature references update");

    my $command = <<PLSQL;
BEGIN
  RNC_UPDATE.update_literature_references();
END;
PLSQL

    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("Literature references update failed " . $DBI::errstr);

    $self->log_plsql_output();
    $self->{'logger'}->info("Literature references update complete");
}


=head2 run_pl_sql_update

    Update Oracle database, import sequences, create new xrefs.

=cut

sub run_pl_sql_update {
    my $self = shift;

    $self->{'logger'}->info("Launching PL/SQL update");

    if ( $self->{'release_type'} ne 'F' and
         $self->{'release_type'} ne 'I' ) {
        $self->{'logger'}->logdie('Incorrect release type');
    }

    my $command = <<PLSQL;
BEGIN
  RNC_UPDATE.NEW_UPDATE('$self->{'release_type'}');
END;
PLSQL

    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("PL/SQL update failed " . $DBI::errstr);

    $self->log_plsql_output();
    $self->{'logger'}->info("PL/SQL update complete");
}


=head2 log_plsql_output

    Auxiliary function for logging PL/SQL output.

=cut

sub log_plsql_output {
    my $self = shift;
    my $output = $self->{'dbh'}->func( 'dbms_output_get' );
    $self->{'logger'}->info($output);
    print $output, "\n";
}


1;