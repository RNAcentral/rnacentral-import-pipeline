#!/usr/bin/perl

# Script for importing embl files with ncRNAs into a database.

use warnings;
use strict;

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::SeqStats;
use Getopt::Long;
use Config::Simple;

use SWISS::CRC64;
use Digest::MD5 qw(md5_hex);
use DBI;
use DBD::Oracle;

# initialize parameters from the config file
my $cfg = new Config::Simple('config.cfg');

my $host    = $cfg->param("db.host");
my $user    = $cfg->param("db.user");
my $pass    = $cfg->param("db.pass");
my $port    = $cfg->param("db.port");
my $sid     = $cfg->param("db.sid");
my $table   = $cfg->param("db.table");
my $dir     = $cfg->param("general.dir");
my $verbose = $cfg->param("general.verbose");

# get command line options to override the config file
&GetOptions (
        'dir:s'  => \$dir,
        'vb:s'   => \$verbose,
        );

&Usage if ($dir eq '');

# locate input files
my @files = get_files($dir);

# connect to the database
my $dsn = "dbi:Oracle:host=$host;sid=$sid;port=$port";
my $dbh = DBI->connect($dsn, $user, $pass) or die( $DBI::errstr . "\n" );

# prepare the table
my $command = "TRUNCATE TABLE $table";
$dbh->do($command);
my $sth = $dbh->prepare("INSERT INTO $table (crc64, len, seq_short, seq_long, ac, version, taxid, md5) VALUES (?, ?, ?, ?, ?, ?, ?, ?)")
             or die "prepare failed: " . $dbh->errstr();  

# loop over files
foreach my $file (@files) {
	
	# read file using Bioperl
	my $stream = Bio::SeqIO->new(-file => $file, -format => 'EMBL');

	my $i = 0;
	while ( (my $seq = $stream->next_seq()) ) {

		$i++;
		print "$i, $seq->display_id\n";		

		my $md5 = md5_hex($seq->seq);	
		my $crc64 = SWISS::CRC64::crc64($seq->seq);

		if ($seq->length < 4000) {
			$sth->execute($crc64, $seq->length, $seq->seq, 'NULL', $seq->display_id, 1, 1, $md5) or die "execute failed: " . $sth->errstr(); 
		} else {
			$sth->execute($crc64, $seq->length, 'NULL', $seq->seq, $seq->display_id, 1, 1, $md5) or die "execute failed: " . $sth->errstr(); 
		}

		if ($verbose eq 'y') {
			print 'seq_version: SV' . $seq->seq_version  . "\n";
			if ( $seq->is_circular ) {
				print "circular\n";
			} else {
				print "linear\n";
			};	
			print 'molecule: ' . $seq->molecule . "\n";
			print 'namespace: ' . $seq->namespace;
			print 'division: ' . $seq->division . "\n";	
			print 'length: BP' . $seq->length . "\n";	
			print 'primary_id: ' . $seq->primary_id . "\n";
			print 'get_dates: ' . join(", ", $seq->get_dates) . "\n";
			print 'desc: ' . $seq->desc . "\n";
			print 'accession_number: ' . $seq->accession_number . "\n";
			print 'get_secondary_accessions: ' . $seq->get_secondary_accessions . "\n";
			print 'keywords: ' . $seq->keywords . "\n";
			print 'species binomial: ' . $seq->species->binomial() . "\n";
			print 'organelle: ' . $seq->species->organelle . "\n";
			print 'common_name: ' . $seq->species->common_name . "\n";
			print 'classification: ' . join(', ', $seq->species->classification) . "\n";					

			print 'seq: ' . $seq->seq . "\n";	
			my $hash_ref = Bio::Tools::SeqStats->count_monomers($seq);
			foreach my $base (sort keys %$hash_ref) {
				print "Number of bases of type ", $base, "= ", %$hash_ref->{$base},"\n";
			}			

			print "\nFeature table:\n";
			for my $feat_object ($seq->get_SeqFeatures) {          

				print 'Location: ' . $feat_object->location->to_FTstring() . "\n"; 
				print 'Start: ' . $feat_object->location->start . "\n"; 
				print 'Stop: ' . $feat_object->location->end . "\n"; 		

				print "primary tag: ", $feat_object->primary_tag, "\n";          
				for my $tag ($feat_object->get_all_tags) {             
					print "  tag: ", $tag, "\n";             
					for my $value ($feat_object->get_tag_values($tag)) {                
						print "    value: ", $value, "\n";             
					}          
				}       
			}

			print "\nAnnotations:\n";
			
			my $anno_collection = $seq->annotation;
			for my $key ( $anno_collection->get_all_annotation_keys ) {
			   my @annotations = $anno_collection->get_Annotations($key);
			   for my $value ( @annotations ) {
			      print "tagname : ", $value->tagname, "\n";
			      # $value is an Bio::Annotation, and also has an "as_text" method
			      print "  annotation value: ", $value->display_text, "\n";	      
			   }
			}	

		} # end if verbose
	} # end while seq
} # end for file

$sth->finish();
$dbh->disconnect();


### SUBS ###

sub get_files {
	opendir(DIR, shift @_) or die $!;

	my @files = ();
	while (my $file = readdir(DIR)) {

	    # ignore files beginning with a period
	    next if ($file =~ m/^\./);
	    next if ($file !~ m/\.ncRNA$/);

		push @files, $dir . '/' . $file;
	}

	closedir(DIR);

	print 'Found ' . scalar @files . " files\n";

	return @files;
}

sub Usage {
    print <<EOF
 $0 --dir <path to input files> [--vb <y/n>]
        
        Mandatory
            --dir     Location of input files 
                        
        Optional       
            --vb      Verbose mode
EOF
;
    exit(1);
}
