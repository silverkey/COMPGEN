#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

my $USAGE = "\nperl $0 [fasta file] [dbname] [create 'Y' or 'N']\n\n";

die $USAGE unless -e $ARGV[0];
die $USAGE unless $ARGV[1];
die $USAGE unless $ARGV[2];

my $fa_file = $ARGV[0];
my $dbname = $ARGV[1];
my $create = 0;
$create = 1 if $ARGV[2] eq 'Y';

my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                         -user => 'mysql_dev',
                                         -pass => 'riiGbs',
                                         -dsn => "dbi:mysql:$dbname",
                                         -create => $create);

my $seqio = Bio::SeqIO->new(-format => 'fasta',
                            -file => $fa_file);

while(my $seq = $seqio->next_seq) {
  $db->insert_sequence($seq->id,$seq->seq);
}
