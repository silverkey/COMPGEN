#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

#use lib '/Users/Remo/src/BioPerl-1.6.1';
#use lib '/Users/Remo/src/autocne';

use Bio::DB::SeqFeature::Store;

my $USAGE = "\nperl $0 [sqlite file] [directory]\n\n";
die $USAGE unless -d $ARGV[1];
chdir($ARGV[1]);
die $USAGE unless -e $ARGV[0];
die $USAGE unless $ARGV[0] =~ /\.sqlite$/;

my $sql_file = $ARGV[0];

my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::SQLite',
                                         -dsn => $sql_file);

my @features = $db->get_features_by_type('match');

foreach my $feature(@features) {

  my $identity = $feature->score;
	
  $identity =~ s/^(\d+\.\d)\d+$/$1/;
	
  my $length = $feature->length;
	
  next unless $length >= 150 && $identity >= 100;
	
  #my @feats = $db->get_features_by_location($feature->seq_id,$feature->start,$feature->end);
  #my $c = 0;
  #foreach my $f(@feats) {
  #  $c ++ if $f->type =~ /^exon/;
  #  $c ++ if $f->type =~ /^match/;
  #}

  #my $seq = $db->fetch_sequence($feature->seq_id,$feature->start,$feature->end);

  #print $feature->seq_id.':'.$feature->start.'-'.$feature->end." $seq\n" unless $c;

  print $feature->seq_id.':'.$feature->start.'-'.$feature->end."\n";
}


#sub get_second_organism_mapping {	
#my $f = shift;
	
