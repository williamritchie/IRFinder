#!/usr/bin/env perl
use strict;
#use Fcntl;

#fcntl(stdout, F_SETPIPE_SZ, 1048576);
#fcntl(fileno(STDOUT), 1031, 1048576);

my $generatedCount = 0;

my $readLen = $ARGV[0];
my $stride = $ARGV[1];

my $lastDirection = 0;

sub reverse_complement {
        my $dna = shift;

        # reverse the DNA sequence
        my $revcomp = reverse($dna);

        # complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}  


my @error = (
{'A' => 'G', 'G' => 'T', 'T' => 'C', 'C' => 'A', 'N' => 'N' },
{'A' => 'T', 'G' => 'C', 'T' => 'A', 'C' => 'G', 'N' => 'N' },
{'A' => 'C', 'G' => 'A', 'T' => 'G', 'C' => 'T', 'N' => 'N' }
);

my $readCount = 0;

sub processRead( $ $ $ ) {
  my $read = shift;
  my $pos = shift;
  my $chr = shift;

  $readCount++;

  my $numN = $read =~ tr/N/N/;
  if ($numN * 2 < $readLen) {
    #only output reads where less than half of the read will be NNNNN

    # generate a single base error in a deterministic manner.
    substr($read,35,1) = $error[$readCount % 3]{substr($read,35,1)};

    if ($lastDirection == 0) {
      print ">RF!$chr!$pos\n";
      print "$read\n";
      $lastDirection = 1;
    }else{
      print ">RR!$chr!$pos\n";
      print reverse_complement($read) . "\n";
      $lastDirection = 0;
    }

  }
}

sub processBuffer( $ $ $ ) {
  my $b = shift;
  my $pos = shift;
  my $chr = shift;

  #while (length($$b) >= $readLen + $stride) {
  while (length($$b) >= $readLen) {
    processRead(substr($$b,0,$readLen), $pos, $chr);
    #my $thisread = substr($$b,0,$readLen);
    $$b = substr($$b,$stride);
    $pos = $pos + $stride;
  }
  return $pos;
}

my $count = 0;
my $chr = '';
my $pos = 1;
my $buffer = '';

while(<STDIN>) {
  chomp;
  $count ++;
  if (m/^>/) {
    s/ .*$//;
    s/^>//;
    $chr = $_;
    $pos = 1;
    $buffer = '';
  }
  else{
    # Should allow into the buffer only valid letters.
    $_ = uc($_);
    s/[^ATCGN]/N/g;
    $buffer .= $_;
    $pos = processBuffer(\$buffer, $pos, $chr);
  }
#  if ($count > 10000) { exit; }
}
