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

    Bio::RNAcentral::SqlldrImportBase

=head1 DESCRIPTION


=cut

package Bio::RNAcentral::SqlldrImportBase;

use strict;
use warnings;

use File::Spec;

use base ('Bio::RNAcentral::OracleUpdate');


sub new {
    my ($class, $opt, $prefix) = @_;

    # run parent constructor
    my $self = $class->SUPER::new($opt);

    my $path;

    if ( $prefix eq 'ac_info' )  {
        $path = $self->get_info_path()
    } elsif ( $prefix eq 'composite_ids' ) {
        $path = $self->get_comp_id_path();
    } elsif ( $prefix eq 'refs' ) {
        $path = $self->get_refs_path();
    } elsif ( $prefix eq 'genome_coordinates' ) {
        $path = $self->get_genomic_locations_path();
    } else {
        $self->{'logger'}->logdie('Incorrect prefix parameter')
    }

    $self->{'local'} = {
        path    => $path,
        ctlfile => File::Spec->catfile($path, "$prefix.ctl"),
        allfile => File::Spec->catfile($path, "$prefix.dat"),
        badfile => File::Spec->catfile($path, "$prefix.bad"),
        logfile => File::Spec->catfile($path, "$prefix.log"),
    };

    return $self;
}


=head2 load_all_references

    Main subroutine for loading the data.

=cut

sub load_all {
}


=head2 _make_ctl_file

    Create a control file used by sqlldr.

=cut

sub _make_ctl_file {
}


=head2 _get_sqlldr_command

    Construct the sqlldr system command.

=cut

sub _get_sqlldr_command {
    my $self = shift;

    my $command = 'sqlldr ' .
                   $self->{'user'} . '/' . $self->{'password'} . '@' .
                  '\"\(DESCRIPTION=\(ADDRESS=\(PROTOCOL=TCP\)' .
                  '\(HOST=' . $self->{'host'} . '\)' .
                  '\(PORT=' . $self->{'port'} . '\)\)' .
                  '\(CONNECT_DATA\=\(SERVICE_NAME=' . $self->{'sid'} . '\)\)\)\" ' .
                  'control=' . $self->{'local'}{'ctlfile'} . ' ' .
                  'bad='     . $self->{'local'}{'badfile'} . ' ' .
                  'log='     . $self->{'local'}{'logfile'} . ' ' .
                  'direct=true errors=99999999 data=' . $self->{'local'}{'allfile'};
    return $command;
}


=head2 _run_sqlldr

    Launch sqlldr and make sure that it runs successfully, log any errors.

=cut

sub _run_sqlldr {
    (my $self, my $command) = @_;
    $self->{'logger'}->info('Launching sqlldr');

    my $status = system($command); # 0 on success
    unless ( $status == 0 ) {
        $self->{'logger'}->logwarn("Couldn't launch sqlldr\n. Command: $command\n Error: $!\n");
    }

    return $status;
}


=head2 _delete_old_log_files

    Delete any .bad and .log files left from previous runs.
    If they are not deleted, it may look as if the current run had an error and not the old one.

=cut

sub _delete_old_log_files {
    my $self = shift;
    my @to_delete = ($self->{'local'}{'badfile'}, $self->{'local'}{'logfile'});
    unlink @to_delete;
}


=head2 _errors_found

    Warn if bad file exists.
    Bad file contains the entries rejected by the database.

=cut

sub _errors_found {
    my $self = shift;

    # TODO: find error messages in the log file

    if (-e $self->{'local'}{'badfile'}) {
        $self->{'logger'}->logwarn("sqlldr import had errors, check $self->{'local'}{'badfile'}");
        return 1;
    } else {
        $self->{'logger'}->info("No bad file");
        return 0;
    }
}


1;
