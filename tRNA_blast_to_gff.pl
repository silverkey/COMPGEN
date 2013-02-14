#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::Tools::GFF;
use Bio::SearchIO;

my $USAGE = "\nperl $0 [blast output file in format 0] [source] [minimum length]\n\n";
die $USAGE unless -e $ARGV[0];
die $USAGE unless $ARGV[1];
die $USAGE unless $ARGV[2];

my $blast_file = $ARGV[0];
my $gff_file = "$blast_file";
$gff_file .= '.gff';

open(OUT,">$gff_file");

my $in = Bio::SearchIO->new(-file => $blast_file,
                            -format => 'blast',
                            -report_format => 0);

while(my $res = $in->next_result) {
  while(my $hit = $res->next_hit) {
    while(my $hsp = $hit->next_hsp) {

      my $f1 = $hsp->feature1;
      my $f2 = $hsp->feature2;

      my $seq_id = $f1->seq_id;

      my $source = $ARGV[1];

      my $type = 'tRNA';

      my $start = $f1->start;

      my $end = $f1->end;

      my $score = $hsp->percent_identity;
      $score =~ s/^(\d+\.\d)\d+$/$1/;

      my $strand = $f2->strand;
      $strand = '+' if $strand eq '1';
      $strand = '-' if $strand eq '-1';

      my $phase = '0';

      print OUT $seq_id."\t".$source."\t".$type."\t".$start."\t".$end."\t".$score."\t".$strand."\t".$phase."\n" if $f1->length >= $ARGV[2];
    }
  }
}

close(OUT);
