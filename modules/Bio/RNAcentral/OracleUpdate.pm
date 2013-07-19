#!/usr/bin/env perl

package Bio::RNAcentral::OracleUpdate;

=pod

=head1 NAME



=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

use strict;
use warnings;

use DBI;
use DBD::Oracle;

use base ('Bio::RNAcentral::Base');


sub new {
    my ($class, $opt) = @_;

    # run parent constructor
    my $self = $class->SUPER::new();

    if ( defined($opt->{'user'}) ) {
        $self->{'opt'}{'user'} = $opt->{'user'};
    }
    if ( defined($opt->{'password'}) ) {
        $self->{'opt'}{'password'} = $opt->{'password'};
    }
    if ( defined($opt->{'host'}) ) {
        $self->{'opt'}{'host'} = $opt->{'host'};
    }
    if ( defined($opt->{'sid'}) ) {
        $self->{'opt'}{'sid'} = $opt->{'sid'};
    }
    if ( defined($opt->{'port'}) ) {
        $self->{'opt'}{'port'} = $opt->{'port'};
    }

    $self->db_oracle_connect();

    return $self;
}


=head2 db_oracle_connect

    Establish a database connection and keep it open.

=cut

sub db_oracle_connect {
    my $self = shift;

    my $dsn = "dbi:Oracle:host=$self->{'opt'}{'host'};sid=$self->{'opt'}{'sid'};port=$self->{'opt'}{'port'}";
    my $dbh = DBI->connect($dsn, $self->{'opt'}{'user'}, $self->{'opt'}{'password'})
              or $self->{'logger'}->logdie( $DBI::errstr . "\n" );

    $self->{'dbh'} = $dbh;
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


=head2 truncate_staging_table

    The staging table should be empty at the beginning of the import.

=cut

sub truncate_staging_table {
    my $self = shift;

    my $command = "TRUNCATE TABLE $self->{'opt'}{'staging_table'}";
    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("Couldn't truncate the staging table");

    $self->{'logger'}->info("Staging table truncated");
}


=head2 create_new_release

    TODO: figure out how many staging tables are going to be used.
    TODO: replace placeholder parameters with real dates etc.

=cut

sub create_new_release {
    my $self = shift;

    my $command = <<SQL;
INSERT INTO $self->{'opt'}{'release_table'}
VALUES(
  (SELECT count(*)+1 FROM rnc_release),
  1,
  '12-JUL-13',
  'F',
  'L',
  '12-JUL-13',
  'perl',
  'test',
  'N'
)
SQL

    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("Couldn't create new release");

    $self->{'logger'}->info("New release created");
}


=head2 run_pl_sql_update

    Update Oracle database, import sequences, create new xrefs.

=cut

sub run_pl_sql_update {
    my $self = shift;

    $self->{'logger'}->info("Launching PL/SQL update");

    my $command = <<PLSQL;
BEGIN
  RNC_UPDATE.NEW_UPDATE();
END;
PLSQL

    $self->{'dbh'}->do($command)
        or $self->{'logger'}->logdie("PL/SQL update failed");

    $self->{'logger'}->info("PL/SQL update complete");
}


1;