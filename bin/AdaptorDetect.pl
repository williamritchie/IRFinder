#!/usr/bin/perl
use strict;
use List::Util qw(sum max min);
use Data::Dumper;

my $DEBUG = 0;
my $maxReadsDetection = 250000;
my $minInsertOutput = 30;

if (scalar @ARGV != 2) {
  print STDERR "Usage: AdaptorDetect.pl reads_1.fastq[.gz] reads_2.fastq[.gz]\n";
  exit(1);
}

my ($in1, $in2) = @ARGV;
my $foundadapt1 = "";
my $foundadapt2 = "";

print "Adaptor Detect v0.7\n";

sub opencmd($$) {
  my ($file, $rw) = @_;
  if ($file =~ m/\.gz/) {
    $file = "gzip -dc < '$file' |" if $rw eq 'r';
    $file = "| gzip > '$file' " if $rw eq 'w';
  }elsif ($file =~ m/\.bz2?/) {
    $file = "bzip2 -dc < '$file' |" if $rw eq 'r';
    $file = "| bzip2 > '$file' " if $rw eq 'w';  
  }elsif ($file =~ m/\.lzop/) {
    $file = "lzop -dc < '$file' |" if $rw eq 'r';
    $file = "| lzop > '$file' " if $rw eq 'w';  
  }else{
    $file = "<$file" if $rw eq 'r';
    $file = ">$file" if $rw eq 'w';
  }
  return $file;
}


open (F1, opencmd($in1, 'r'));
open (F2, opencmd($in2, 'r'));

my $mininsert = 20;

my ($L1, $L2);

my $readlen = 0;
my $consideredReads = 0;


my (%A1_10, %A1_10_frag_dist, %A1_10_score, %A2_10, %A2_10_frag_dist, %A2_10_score);

my $readcount = 0;
while (<F1>) {
	chomp(my $ID1 = $_);
	chomp(my $ID2 = <F2>);
	chomp(my $Seq1 = <F1>);
	chomp(my $Seq2 = <F2>);
	chomp(my $IDq1 = <F1>);
	chomp(my $IDq2 = <F2>);
	chomp(my $Qual1 = <F1>);
	chomp(my $Qual2 = <F2>);

        $readcount ++;

    $readlen = "[" . length($Seq1) . "," . length($Seq2) . "]";

    $Seq1 =~ s/N*$//i;
    $Seq2 =~ s/N*$//i;    
	next if ($Seq1 =~ m/NNNNNN/i || $Seq2 =~ m/NNNNNN/i);
	next if (length($Seq1) < 30 || length($Seq2) < 30);

	$consideredReads ++;

	my $search1 = $Seq1;
	my $search2 = reverse($Seq2);
	$search1 =~ tr/ATCGatcgNn/01230123../;	
	$search2 =~ tr/ATCGatcgNn/10321032++/;	


  my $best_score = -10000;
  my @scores;
  my $best_insert = 10000;
  my $best_overlap = 10000;
  my $len1 = length($search1);
  my $len2 = length($search2);
  my $maxinsert = max($len1, $len2) - 1;
  for (my $insert = $mininsert; $insert <= $maxinsert; $insert++) {
    my $xor = "";
    my $overlap = $insert;
    if ($len2 == $len1) {
      $xor = substr($search1, 0, $insert) ^ substr($search2, -$insert);  
    }elsif ($len2 > $len1) {
      $overlap = min($insert, $len1);
      $xor = substr($search1, 0, $insert) ^ substr($search2, -$insert, $overlap);      
    }else{
      $overlap = min($insert, $len2);
      $xor = substr($search1, $insert-$overlap, $overlap) ^ substr($search2, -$overlap);          
    }
    
    #my $xor = substr($search1, $insert-$overlap, $overlap) ^ reverse(substr($search2, 0, $i));

    my $mismatch = $xor =~ tr/\1\2\3/XXX/; ## 0/1 1/2 1/3 .. 2/3
    my $match = $xor =~ tr/\0//;  ## 0/0 1/1 2/2 3/3
    ## Note no match for ./. because . has been translated to + in one sample.
    my $score = ($match - $mismatch)/$overlap;

    push @scores, $score;
    #if ($score > 0.3 && $mismatch < ($overlap / 6)) {
        if ($score > $best_score) {
          $best_score = $score;
          $best_insert = $insert;
          $best_overlap = $overlap;
        }
    #}
  }

  my $frag_len = 1000;
  if ($best_insert != 10000) {
   my $second_best_score = 0.3;
   $second_best_score = (sort @scores)[-2] if ((scalar @scores) >= 2);
   if ($best_score > 0.7 && $best_score >= ($second_best_score+0.4)) {

    my $adapt1 = "";
    my $adapt2 = "";
    $adapt1 = substr($Seq1,$best_insert) if ($len1 > $best_insert);
    $adapt2 = substr($Seq2,$best_insert) if ($len2 > $best_insert);
    
    if (! ($adapt1 =~ m/N/i) ) {
      $A1_10{substr($adapt1,0,10)}++ if (length($adapt1)>=10);
      $A1_10_frag_dist{substr($adapt1,0,10)}{$best_insert}++ if (length($adapt1)>=10);
    }
    if (! ($adapt2 =~ m/N/i) ) {
      $A2_10{substr($adapt2,0,10)}++ if (length($adapt2)>=10);
      $A2_10_frag_dist{substr($adapt2,0,10)}{$best_insert}++ if (length($adapt2)>=10);
    }

    if ($DEBUG) {
if (length($adapt1)>=6 && substr($adapt1,0,6) ne "AGATCG") {
      print sprintf("%-4s",($len1-$best_insert)) . " " . sprintf("%0.2f",$best_score) . "  " . sprintf("%-75s",substr($Seq1,$best_insert)) . substr($Seq2,$best_insert) . "\n"; 
      print( (" " x ($len2-$best_insert)) . $Seq1 . "\n");
      my $Seq2Comp = $Seq2;
      $Seq2Comp =~ tr/ATCGatcg/TAGCtagc/;
      print reverse($Seq2Comp) . "\n";
      print reverse($Seq2) . "\n";
}
    }
   }
  }

  if ($readcount >= $maxReadsDetection) {
    last;
  }
}

close F1;
close F2;

for my $adapt (keys %A1_10_frag_dist) {
  my $score = 0;
  for my $insert (keys %{$A1_10_frag_dist{$adapt}}) {
    $score += log(1+$A1_10_frag_dist{$adapt}{$insert})/log(1.5);
  }
  $A1_10_score{$adapt} = $score;
}
for my $adapt (keys %A2_10_frag_dist) {
  my $score = 0;
  for my $insert (keys %{$A2_10_frag_dist{$adapt}}) {
    $score += log(1+$A2_10_frag_dist{$adapt}{$insert})/log(1.5);
  }
  $A2_10_score{$adapt} = $score;
}

print "--Detection detail--\n";

my $tot = 0;
my $c = 0;
my $keyone = "";
my $one = 0;
my $next = 0;
for my $key ( sort {$A1_10_score{$b}<=>$A1_10_score{$a}} keys %A1_10_score) {
           print (sprintf("%10.3f  ", $A1_10_score{$key} ), sprintf("%5s  ", scalar keys %{$A1_10_frag_dist{$key}} ) , sprintf("%5s  ", ($A1_10{$key})), $key, "\n")  if (($c < 4) && ($A1_10{$key} > 1));
           $keyone = $key if ($c == 0);
           $one = $A1_10_score{$key} if ($c == 0);
           $next = $A1_10_score{$key} if ($c == 1);
		   $tot += $A1_10{$key};
		   $c ++;
}
#print "Keys: " . keys(%A1_10) . "  Total: " . $tot . "\n";
if ( (($next * 1.5)<$one) && ($one > 8)) {
  print "FOUND-A1 - $readlen - $one / $next:  $keyone\n";
  $foundadapt1 = $keyone;
  print "Adaptor1: $foundadapt1\n";  
}else{
  print "MISS-A1 - $readlen - $one / $next:  $keyone\n"; 
}

$tot = 0;
$c = 0;
$keyone = "";
$one = 0;
$next = 0;
for my $key ( sort {$A2_10_score{$b}<=>$A2_10_score{$a}} keys %A2_10_score) {
           print (sprintf("%10.3f  ", $A2_10_score{$key} ), sprintf("%5s  ", scalar keys %{$A2_10_frag_dist{$key}} ) , sprintf("%5s  ", ($A2_10{$key})), $key, "\n")  if (($c < 4) && ($A2_10{$key} > 1));
           #print (sprintf("%5s  ", ($A2_10{$key})), $key, "\n")  if (($c < 5) && ($A2_10{$key} > 1));
           $keyone = $key if ($c == 0);
           $one = $A2_10_score{$key} if ($c == 0);
           $next = $A2_10_score{$key} if ($c == 1);
		   $tot += $A2_10{$key};
		   $c ++;
}
#print "Keys: " . keys(%A2_10) . "  Total: " . $tot . "\n";
if ( (($next * 1.5)<$one) && ($one > 8)) {
  print "FOUND-A2 - $readlen - $one / $next:  $keyone\n";
  $foundadapt2 = $keyone;
  print "Adaptor2: $foundadapt2\n";  
}else{
  print "MISS-A2 - $readlen - $one / $next:  $keyone\n"; 
}

print "Total reads considered during auto-detect: $consideredReads\n";




sub doTrim($$$$) {
    my ($Seq1, $Seq2, $a1, $a2) = @_;

    my $returnInsert = -1;

    $Seq1 =~ s/N*$//i;
    $Seq2 =~ s/N*$//i;

	my $search1 = $a1 . $Seq1;
	my $search2 = reverse($a2 . $Seq2);
	$search1 =~ tr/ATCGatcgNn/01230123../;
	$search2 =~ tr/ATCGatcgNn/10321032++/;

  my $best_score = -10000;
  my @scores;
  my $best_insert = 10000;
  my $lR1 = length($Seq1);
  my $lR2 = length($Seq2);
  my $lA1 = length($a1);
  my $lA2 = length($a2);
  
  my $maxinsert = max($lR1, $lR2) - 1;
  for (my $insert = 0; $insert <= $maxinsert; $insert++) {
    my $overlap = min($lR1, $insert+$lA2) - max(-$lA1, $insert-$lR2);
    my $xor = substr($search1, max($lA1+$insert-$lR2,0), $overlap) ^ substr($search2, max($lR2-$insert-$lA1, 0), $overlap);
    
    my $match = $xor =~ tr/\0/X/;  ## 0/0 1/1 2/2 3/3
    my $mismatch = $xor =~ tr/\1\2\3//; ## 0/1 1/2 1/3 .. 2/3
    ## Note no match for ./. because . has been translated to + in one sample.
    my $score = ($match - $mismatch)/$overlap;

    push @scores, $score;
        if ($score > $best_score) {
          $best_score = $score;
          $best_insert = $insert;
          #$best_overlap = $overlap;
        }
  }

  if ($best_insert != 10000) {
   my $second_best_score = 0.3;
   $second_best_score = (sort @scores)[-2] if ((scalar @scores) >= 2);
   if ($best_score > 0.7 && $best_score >= ($second_best_score+0.4)) {
    $returnInsert = $best_insert;

   }
  }
  return $returnInsert;
}


print "===SUMMARY===\n";
if ((length($foundadapt1) + length($foundadapt2)) < 10) {
  print "Fatal: Neither adaptor positively detected\n";
  print STDERR "Fatal: Neither adaptor positively detected\n";
  exit 1;
}
if ((length($foundadapt1) + length($foundadapt2)) < 20) {
  print "Fatal: Only one adaptor was positively detected\n";
  print STDERR "Fatal: Only one adaptor was positively detected\n";
  exit 1;
}
if ($foundadapt1 ne 'AGATCGGAAG') {
  print "Warning: Adaptor 1 is not the very common Illumina TruSeq Adaptor. Is this a specialised experiment? Do you expect a non-standard adaptor?\n";
  print STDERR "Warning: Adaptor 1 is not the very common Illumina TruSeq Adaptor. Is this a specialised experiment? Do you expect a non-standard adaptor?\n";
  exit 1;
}
if ($foundadapt2 ne 'AGATCGGAAG') {
  print "Warning: Adaptor 2 is not the very common Illumina TruSeq Adaptor. Is this a specialised experiment? Do you expect a non-standard adaptor?\n";
  print STDERR "Warning: Adaptor 2 is not the very common Illumina TruSeq Adaptor. Is this a specialised experiment? Do you expect a non-standard adaptor?\n";
  exit 1;
}
print "Standard Illumina TruSeq Adaptor 1 and 2 detected: $foundadapt1, $foundadapt2\n";
