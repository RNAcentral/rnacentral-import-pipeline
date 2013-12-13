#!/usr/bin/env perl

package Bio::RNAcentral::UpdateRefsUsingEuropePMC;

=pod

=head1 NAME

    Bio::RNAcentral::UpdateRefsUsingEuropePMC

=head1 DESCRIPTION



=cut

use strict;
use warnings;

use LWP::Simple qw(get);
use JSON        qw(from_json);
use URI::Escape;


use base ('Bio::RNAcentral::OracleUpdate');


sub new {
    my ($class, $opt, $path, $prefix) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt);

    return $self;
}


sub update_entries {
    my $self = shift;

    my $sth = $self->{'dbh'}->prepare('SELECT id, title FROM rnc_references WHERE pmid is null') or die "Couldn't prepare statement: " . $self->{'dbh'}->errstr;

    $sth->execute() or die "Couldn't execute statement: " . $sth->errstr;

    $self->{'update_sth'} = $self->{'dbh'}->prepare('UPDATE rnc_references SET pmid = ?, pmcid = ?, epmcid = ?, doi = ? WHERE id = ?')
      or die "Couldn't prepare statement: " . $self->{'dbh'}->errstr;

    while (my $data = $sth->fetchrow_hashref()) {
        # use Data::Dumper;
        # print Dumper($data);

        if ( defined($data->{'TITLE'}) ) {
            print "$data->{'TITLE'}\n";
            my $epmc_data = $self->_get_europepmc_data($data->{'TITLE'});
        use Data::Dumper;
        print Dumper($epmc_data);
            if ($epmc_data) {
                $self->_save_europepmc_data($epmc_data, $data->{'ID'});
            }
        }
    }
}



sub _get_europepmc_data {

    my ($self, $title) = @_;

    my $query = "query=TITLE:" . uri_escape($title) . "&format=json&resulttype=lite";
    my $url   = "http://www.ebi.ac.uk/europepmc/webservices/rest/search/" . $query;

    my ($response, $data);

    $response = get($url);

    if (defined($response)) {
        $data = from_json($response);
    } else {
        print "Problem fetching $url";
    }

    if ( $data->{'hitCount'} > 0 ) {
        return {
            'epmcid' => $self->_get_attribute($data, 'id'),
            'pmid '  => $self->_get_attribute($data, 'pmid'),
            'pmcid'  => $self->_get_attribute($data, 'pmcid'),
            'doi'    => $self->_get_attribute($data, 'doi'),
        }
    } else {
        print "No match for $title\n";
        return ();
    }
}


sub _save_europepmc_data {

    my ($self, $data, $row_id) = @_;

    my @array = ($data->{'pmid'}, $data->{'pmcid'}, $data->{'epmcid'}, $data->{'doi'}, $row_id);
    $self->{'update_sth'}->execute_array(@array)
      or die "Couldn't execute statement: " . $self->{'update_sth'}->errstr;
}


sub _get_attribute {
    my ($self, $data, $attribute) = @_;

    if ( exists $data->{'resultList'}{'result'}[0]{$attribute} ) {
        return $data->{'resultList'}{'result'}[0]{$attribute};
    } else {
        return undef;
    }
}


1;
