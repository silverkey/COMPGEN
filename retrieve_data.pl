#!/usr/bin/perl
use strict;
use warnings;
use LWP::Simple;

my $query = 'biomol_rRNA[PROP]';
#assemble the esearch URL
my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=nucleotide&term=$query&usehistory=y";
#post the esearch URL
my $output = get($url);

#parse WebEnv, QueryKey and Count (# records retrieved)
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

#open output file for writing
open(OUT, ">rRNA.fa") || die "Can't open file!\n";

print "Downloading $count sequences....\n\n";

#retrieve data in batches of 500
my $retmax = 5000;
for (my $retstart = 0; $retstart <= $count; $retstart += $retmax) {
  my $efetch_url = $base ."efetch.fcgi?db=nucleotide&WebEnv=$web";
  $efetch_url .= "&query_key=$key&retstart=$retstart";
  $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
  my $efetch_out = get($efetch_url);
  while(! defined($efetch_out)){
    sleep 1;
    $efetch_out = get($efetch_url);
  }
  print "printing $retstart....\n";
  print OUT "$efetch_out";
  sleep 3;
}
close OUT;

print "THE END\n";
