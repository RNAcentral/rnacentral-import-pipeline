#!/usr/bin/env perl

package Bio::RNAcentral::UpdateRefsUsingEuropePMC;

=pod

=head1 NAME

    Bio::RNAcentral::UpdateRefsUsingEuropePMC

=head1 DESCRIPTION

    Import missing bibliographic information from EuropePMC.
    To re-update the data:
    UPDATE rnc_references SET epmc_updated = 0;

=cut

use strict;
use warnings;

use LWP::Simple qw(get);
use JSON        qw(from_json);
use URI::Escape;


use base ('Bio::RNAcentral::OracleUpdate');


sub new {
    my ($class, $opt) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt);

    return $self;
}

=head2 update_bibliographic_info

    Main function for importing missing data from EuropePMC.

=cut

sub update_bibliographic_info {
    my $self = shift;

    $self->{'logger'}->info('Updating bibliographic info');

    $self->_prepare_update_statements();
    my $missing_data = $self->_get_entries_with_missing_data();

    my $updated = my $total = 0;
    while (my $data = $missing_data->fetchrow_hashref()) {
        $updated += $self->_update_entry($data);
        $total += 1;
    }
    $self->{'logger'}->info("$total entries with missing data, $updated updated entries.");
    $self->{'logger'}->info("Bibliographic info update completed.");
}


=head2 _prepare_update_statements

    Preparing the update statement improves performance
    since the stored statement can be reused multiple times.

=cut

sub _prepare_update_statements {

    my $self = shift;
    my $cmd = 'UPDATE rnc_references SET pmid = ?, pmcid = ?, epmcid = ?, doi = ?, epmc_updated = 1 WHERE id = ?';
    $self->{'save_new_data'} = $self->{'dbh'}->prepare_cached($cmd)
        or $self->{'logger'}->logdie("Couldn't prepare statement: " . $self->{'dbh'}->errstr);

    $cmd= 'UPDATE rnc_references SET epmc_updated = 1 WHERE id = ?';
    $self->{'mark_as_done_sth'} = $self->{'dbh'}->prepare_cached($cmd)
        or $self->{'logger'}->logdie("Couldn't prepare statement: " . $self->{'dbh'}->errstr);
}


=head2 _get_entries_with_missing_data

    Get result handle of the table rows containing missing
    bibliographic information.

    Analyze entries without pubmed id or PMC/EuropePMC id or DOI
    which have meaningful titles.

=cut

sub _get_entries_with_missing_data {

    my $self = shift;

    # sql query to retrieve incomplete entries
    my $cmd = <<SQL;
    SELECT id, title, pmid FROM rnc_references
    WHERE (epmc_updated is null or epmc_updated = 0)
    AND (location != 'Unpublished.' AND location NOT LIKE 'Submitted%')
SQL

    my $sth = $self->{'dbh'}->prepare($cmd)
        or $self->{'logger'}->logdie("Couldn't prepare statement: " . $self->{'dbh'}->errstr);
    $sth->execute() or $self->{'logger'}->logwarn("Couldn't execute statement: " . $sth->errstr);

    return $sth;
}


=head2 _update_entry

    Retrieve and store new bibliographic data.

=cut

sub _update_entry {
    my ($self, $missing_data) = @_;
    my $updated = 0;
    if ( defined($missing_data->{'TITLE'}) ) {
        my $epmc_data = $self->_get_europepmc_data($missing_data);

        if ($epmc_data) {
            $self->_save_europepmc_data($epmc_data, $missing_data->{'ID'});
            $updated = 1;
        } else {
            $self->_mark_as_done($missing_data->{'ID'});
        }
    }
    return $updated;
}


=head2 _get_europepmc_data

    Search EuropePMC using REST API.
    Spaces should be replaced by %20 in the title.

=cut

sub _get_europepmc_data {

    my ($self, $missing_data) = @_;
    my ($response, $data, $query);

    $self->{'logger'}->info($missing_data->{'TITLE'});

    # try using pubmed id
    if ( defined($missing_data->{'PMID'}) ) {
        $query = "query=EXT_ID:" . $missing_data->{'PMID'};
    } else { # otherwise search by title
        $query = "query=TITLE:" . uri_escape($missing_data->{'TITLE'});
    }

    $query .= "&format=json&resulttype=lite";
    my $url = "http://www.ebi.ac.uk/europepmc/webservices/rest/search/" . $query;

    $response = get($url);

    if (defined($response)) {
        $data = from_json($response);
    } else {
        $self->{'logger'}->logwarn("Problem fetching $url");
    }

    if ( $data->{'hitCount'} > 0 ) {
        return {
            'epmcid' => $self->_get_result_attribute($data, 'id'),
            'pmid'   => $self->_get_result_attribute($data, 'pmid'),
            'pmcid'  => $self->_get_result_attribute($data, 'pmcid'),
            'doi'    => $self->_get_result_attribute($data, 'doi'),
        }
    } else {
        $self->{'logger'}->info("No match for $missing_data->{'TITLE'}");
        return ();
    }
}


=head2 _get_result_attribute

    Convenience method for accessing the json result object.

=cut

sub _get_result_attribute {
    my ($self, $data, $attribute) = @_;

    if ( exists $data->{'resultList'}{'result'}[0]{$attribute} ) {
        return $data->{'resultList'}{'result'}[0]{$attribute};
    } else {
        return '';
    }
}


=head2 _save_europepmc_data

    Store the newly retrieved data in the database.

=cut

sub _save_europepmc_data {

    my ($self, $data, $row_id) = @_;

    $self->{'logger'}->info("PMID: $data->{'pmid'}, PMCID: $data->{'pmcid'}, EuropePMC: $data->{'epmcid'}, DOI: $data->{'doi'}");

    my @array = ($data->{'pmid'}, $data->{'pmcid'}, $data->{'epmcid'}, $data->{'doi'}, $row_id);
    $self->{'save_new_data'}->execute(@array)
        or $self->{'logger'}->logwarn("Couldn't execute statement: " . $self->{'save_new_data'}->errstr);
}


=head2 _mark_as_done

    Set the "updated" column value to 1 so that this entry
    is not re-analyzed during the next update.

=cut

sub _mark_as_done {
    my ($self, $id) = @_;
    $self->{'mark_as_done_sth'}->execute($id)
        or $self->{'logger'}->logwarn("Couldn't execute statement: " . $self->{'save_new_data'}->errstr);
}


1;
