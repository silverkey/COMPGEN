#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::FeatureIO;

# TO MAKE THIS WORK I HAVE
# MODIFIED A LINE ON BP MODULE
# sudo nano +903 /usr/share/perl5/Bio/FeatureIO/gff.pm

my $USAGE = "\nperl $0 [gff file] [length cutoff]\n\n";

die $USAGE unless -e $ARGV[0];
die $USAGE unless $ARGV[1];

my $gff_file = $ARGV[0];
my $gff_out = "filtered_$ARGV[0]";

my $in = Bio::FeatureIO->new(-file => $ARGV[0],
                             -format => 'gff');

my $out = Bio::FeatureIO->new(-file => ">$gff_out",
                              -format => 'gff');

while(my $f = $in->next_feature) {
  if($f->length >= $ARGV[1]) {
    $out->write_feature($f);
  }
}
