#!/usr/bin/perl

#0	1	2	3						4	5	6	7	8	9	10		11
#1      11868   14409   ENST00000456328/processed_transcript/DDX11L1    0       +       11868   14409   0       3       359,109,1189,   0,744,1352,


open EXON, '>', $ARGV[0];
open INTRON, '>', $ARGV[1];

while (<STDIN>) {
  chomp;
  @f = split /\t/;

  $trans_start = $f[1];

  @length = split /,/, $f[10];
  @start = split /,/, $f[11];
  $chr = $f[0];
  ($gene_id,$gene_name) = $f[3] =~ /\/([^\/]*)\/([^\/]*)$/;
  $dir = $f[5];

  $last_end = undef;
  while (@length) {
    $start = shift @start;
    $length = shift @length;
    if (defined($last_end)) {
        #only output if the intron has length.
        if (($last_end+1) < ($start-1)) {
          print INTRON join("\t", $chr, $trans_start+$last_end, $trans_start+$start, "$gene_name/$gene_id/$dir", $f[4], $f[5]), "\n";
        }
    }
    #print EXON "$chr\t" . ($trans_start+$start) . "\t" . ($trans_start+$start+$length) . "\t$name\n";
    print EXON join("\t", $chr, $trans_start+$start, $trans_start+$start+$length, "$gene_name/$gene_id/$dir", $f[4], $f[5]), "\n";
    $last_end = $start+$length;
  }
}

close INTRON;
close EXON;
