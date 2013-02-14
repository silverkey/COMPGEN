use strict;
use warnings;
use Bio::SeqIO;

my $in = Bio::SeqIO->new(-file => $ARGV[0],
                         -format => 'fasta');

my $out = Bio::SeqIO->new(-file => ">human_$ARGV[0]",
                          -format => 'fasta');


while(my $seq = $in->next_seq) {
  $out->write_seq($seq) if ($seq->id =~ /Homo/);
}
