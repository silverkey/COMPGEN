#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

use Bio::DB::SeqFeature::Store;
use Bio::SeqIO;
use Bio::Seq;

# -----------------------------------------------------------------
# The script is organized in the way to be modulized at some time
# Written on Blast runned with CI as query
#------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# THE KEY OF THE HASH $SPECIES HAVE TO BE IN THE ATTRIBUTE 'O=...' IN THE GFF FROM THE BLAST
# ------------------------------------------------------------------------------------------

my $USAGE = "\nperl $0 [length cutoff] [identity cutoff]\n\n";

die $USAGE unless $ARGV[0];
die $USAGE unless $ARGV[1];

# The species analized
my $species;
$species->{CI} = 'bp_ciona_intestinalis_59';
$species->{CS} = 'bp_ciona_savignyi_59';
$species->{PM} = 'bp_phallusia_1';
$species->{OD} = 'bp_oikopleura_3';

# The maximum number of overlapping hits we can get if we consider one match from each species
my $max_overlapping_aln = (scalar(keys %$species) - 1);

my $organism_string = join('_',keys(%$species));
my $NAME_STRING = "$ARGV[0]\_$ARGV[1]\_$organism_string";
my $DIR_NAME = "aln_$NAME_STRING";
my $FILE_NAME = "selected_$NAME_STRING\.xls";
system("mkdir $DIR_NAME");
chdir("$DIR_NAME");
open(OUT,">$FILE_NAME");

# the species used as query in the blast and which genome db contains all the hits
my $reference = 'CI'; 

# for speed use the genome with the smaller number of blast hits;
my $first = 'OD';
#my $first = 'PM';

# Definition of the databases
my $other_species_db = {};
foreach my $organism(keys %$species) {
  if($organism ne $reference) {
    $other_species_db->{$organism} = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                                                     -dsn => 'dbi:mysql:'.$species->{$organism},
                                                                     -user => 'mysql_dev',
                                                                     -pass => 'riiGbs');
  }
}
my $ref_db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                             -dsn => 'dbi:mysql:'.$species->{$reference},
                                             -user => 'mysql_dev',
                                             -pass => 'riiGbs');

# Count the interesting results
my $nr = 0;

print OUT "Result id";
print OUT "\tRef id\tRef start\tRef end\tStrand\tHit id start end\tHit organism" x scalar(keys %$species);
print OUT "\n";

# Select all the matches between the reference and the first genomes
my @features = $ref_db->features(-type => 'match',
                                 -attributes => {O => $first});

# Filter the matches and search for ovelapping matches in other genomes
foreach my $feature(@features) {

  my @res;

  my @Tattribute = $feature->attributes('T');
  my @Oattribute = $feature->attributes('O');
  my $identity = $feature->score;
  $identity =~ s/^(\d+\.\d)\d+$/$1/;
  my $length = $feature->length;

  # filter on length and identity
  next unless $length >= $ARGV[0] && $identity >= $ARGV[1];

  # Filter out matches overlapping exons
  my @fos1 = $ref_db->features(-seq_id => $feature->seq_id,
                               -start => $feature->start,
                               -end => $feature->end,
                               -type => 'exon');
  next if scalar(@fos1);

  # Filter out matches overlapping tRNAs
  my @fos2 = $ref_db->features(-seq_id => $feature->seq_id,
                               -start => $feature->start,
                               -end => $feature->end,
                               -type => 'tRNA');
  next if scalar(@fos2);

  # What arrives here is a potential result
  push(@res,$feature);

  # Search for overlapping in other species
  foreach my $organism(keys %$species) {
    next if $organism eq $reference;
    next if $organism eq $first;
    my @fs = $ref_db->features(-seq_id => $feature->seq_id,
                               -start => $feature->start,
                               -end => $feature->end,
                               -type => 'match',
                               -attributes => {O => $organism});
    next unless scalar(@fs);

    my $e = 1000000000000;
    my $selected = 0;

    # Take only the match with the highest e-value from each species
    foreach my $f(@fs) {
      next unless $f->length >= $ARGV[0] && $f->score >= $ARGV[1];
      my @e = $f->attributes('E');
      if($e[0] < $e) {
        $selected = $f;
        $e = $e[0];
      }
    }
    next unless $selected;
    push(@res,$selected);
    $selected = 0;
  }
  # Nothing to do if we do not get one match from each species
  next unless scalar(@res) == $max_overlapping_aln;
  
  # -----------------------------------------------
  # We are here means we have a significant result!
  # -----------------------------------------------
  $nr++;
  my $print = "$nr\t";

  # Let's build the longer range on the reference genome
  my $ref_start = 10000000000;
  my $ref_end = 0;
  my $ref_id;

  # The hashref that will contain the seq objects to build the MSA
  my $UCE = {};

  foreach my $res(@res) {
    my @T = $res->attributes('T');
    my @O = $res->attributes('O');
    $print .= $res->seq_id."\t".$res->start."\t".$res->end."\t".$res->strand."\t".$res->score."\t".$T[0]."\t".$O[0]."\t";
    $ref_id = $res->seq_id;

    my @rf = split(/ /,$T[0]);
    my $rdb = $other_species_db->{$O[0]};

    my $sid = "$O[0]\:$rf[0]\:$rf[1]\-$rf[2]\:".$res->strand;
    my $sstring = $rdb->fetch_sequence($rf[0],$rf[1],$rf[2]);
    my $pre_seq = Bio::Seq->new(-id => $sid,
                                -seq => $sstring);
    if($res->strand == -1) {
       $UCE->{$O[0]} = $pre_seq->revcom;
    }
    elsif($res->strand == 1) {
       $UCE->{$O[0]} = $pre_seq;
    }
    else {
      die "\n\nSTRAND ERROR!!!!!\n\n;"
    }

    if($res->start < $ref_start) {
      $ref_start = $res->start;
    }
    if($res->end > $ref_end) {
      $ref_end = $res->end;
    }
  }

  $print =~ s/\t$/\n/;
  print OUT $print;

  my $sid = "$reference\:$ref_id\:$ref_start\-$ref_end\:1";
  my $sstring = $ref_db->fetch_sequence($ref_id,$ref_start,$ref_end);
  $UCE->{$reference} = Bio::Seq->new(-id => $sid,
                                     -seq => $sstring);

  my $fasta = Bio::SeqIO->new(-file => ">$nr.fa",
                              -format => 'fasta');

  foreach my $key(keys %$UCE) {
    my $seq = $UCE->{$key};
    $fasta->write_seq($seq);
  }
}
