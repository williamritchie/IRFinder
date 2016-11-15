#!/usr/bin/perl
use strict;

my %genes;

sub intronNumber
{
  my ($gene) = shift;
  $genes{$gene} ++;
  return $genes{$gene};
}

sub processIntron {
  my ($intron, $overlaps) = @_;
  my ($chr, $start, $end, $gene, $score, $dir) = split /\t/, $intron;

  my $len = $end-$start;
  my $excl = 0;
  my $antisense_dirty = 0;
  my $excluded_by_exon = 0;

  my @intron_seg=();
  push @intron_seg, {'start' => $start, 'end' => $end};

  foreach my $overlap (@$overlaps) {
    # $overlap->{start, end, type}
    if ($overlap->{'type'} eq 'A') {
      $antisense_dirty = 2 if ($antisense_dirty < 2);
      #ignore anti-sense, but mark dirty
    }elsif ($overlap->{'type'} eq 'AE') {
      $antisense_dirty = 1 if ($antisense_dirty < 1);
      #ignore anti-sense, but mark dirty
    }elsif ($overlap->{'type'} eq 'E' && $overlap->{'start'} < $start && $overlap->{'end'} > $end) {
      #print STDERR "Found an exon/feature entirely covering this intron, skpping\n";
      $excluded_by_exon = 1;
    }else{
      # We want to exclude this segment from our intron.

      foreach my $seg (@intron_seg) {
        if ($seg->{'end'}==0) {
          # do nothing, this segment has already been deleted.
        }elsif ($overlap->{'end'} <= $seg->{'start'}) {
          # end is before the start, skip it
        }elsif ($overlap->{'start'} >= $seg->{'end'}) {
          # start is after the end, skip it
        }elsif ($overlap->{'start'} <= $seg->{'start'} && $overlap->{'end'} >= $seg->{'end'}) {
          # exclude entirely covers a segment (equality or beyond) then remove it
          $seg->{'start'}=0; 
          $seg->{'end'}=0; 
        }elsif ($overlap->{'start'} <= $seg->{'start'}) {
          # start is before the start, trim the start
          $seg->{'start'} = $overlap->{'end'};
        }elsif ($overlap->{'end'} >= $seg->{'end'}) {
          # end is after the end, trim the end
          $seg->{'end'} = $overlap->{'start'};
        }else{
          # start inside, end inside - split it
          push @intron_seg, {'start'=>$overlap->{'end'},'end'=>$seg->{'end'}};
          $seg->{'end'}=$overlap->{'start'};
        }
      }
    }
  }
  # Procesed all overlaps.
  # result fragments, in no specific order, in @intron_seg.
  #print $intron, "\n";
  my $newlen = 0;
  my $newstart;
  my $newend = 0;
  my @sizes;
  my @starts;
  foreach my $seg (sort {$a->{'start'} <=> $b->{'start'}} @intron_seg) {
    if ($seg->{'end'} != 0) {
      $newstart = $seg->{'start'} if (!$newstart);
      $newend = $seg->{'end'} if ($seg->{'end'} > $newend);
      push @starts, $seg->{'start'} - $newstart;
      push @sizes, $seg->{'end'}-$seg->{'start'};
    }
    #print join("\t", "", $seg->{'start'}, $seg->{'end'}), "\n";
    $newlen += $seg->{'end'}-$seg->{'start'};
  }
  if ($newlen > 40 && ($newlen/$len) >= 0.7) {
    my $antisense_text = 'clean';
    if ($excluded_by_exon >= 1) {
	$antisense_text = 'known-exon';
        $antisense_text .= '+anti-near' if ($antisense_dirty >= 1);
        $antisense_text .= '+anti-over' if ($antisense_dirty >= 2);
    }else{
        $antisense_text = 'anti-near' if ($antisense_dirty >= 1);
        $antisense_text = 'anti-over' if ($antisense_dirty >= 2);
    }
    print join("\t", $chr, $newstart, $newend, join("/",$gene,intronNumber($gene),$start,$end,$len,$len-$newlen,$antisense_text), $score, $dir, $newstart, $newend, "255,0,0", scalar @sizes, join(",",@sizes), join(",",@starts)), "\n"; 

    if ($len >= 110) {
      print OF50 join("\t", $chr, $start+5, $start+55, "S", 0, $dir, $start+5, $start+55, "255,0,0", 1, 50, 0), "\n";
      print OF50 join("\t", $chr, $end-55, $end-5, "S", 0, $dir, $end-55, $end-5, "255,0,0", 1, 50, 0), "\n";
    }
#    if ($len >= 210) {
#      print OF50 join("\t", $chr, $start+55, $start+105, "E", 0, $dir, $start+55, $start+105, "255,0,0", 1, 50, 0), "\n";
#      print OF50 join("\t", $chr, $end-105, $end-55, "E", 0, $dir, $end-105, $end-55, "255,0,0", 1, 50, 0), "\n";
#    }
    print OF1 join("\t", $chr, $start, $dir), "\n";
    print OF1 join("\t", $chr, $end, $dir), "\n";
  }
}



#### MAIN ####

if (! (scalar @ARGV == 2) ) {
  print STDERR "Usage: cat inputBedIntersection | ./thisTool.pl out2 out3 > out1\n";
  exit(1);
}

open OF50, '>', $ARGV[0];
open OF1, '>', $ARGV[1];

my $lastintron = '';
my @overlaps;
while(<STDIN>) {
  chomp;

## Directional
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       134895  135807  E       0       -       5
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       135135  135900  E       0       -       98
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       135230  136040  X       0       -       238
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       136070  136410  X       0       -       340
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       136440  136710  X       0       -       270
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       136750  137100  X       0       -       350
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       137140  137790  X       0       -       480
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       137615  139384  E       0       -       5
#1       736543  741178  RP11-206L10.8/ENSG00000230092/- 0       -       1       736253  736548  E       0       -       5
#1       736543  741178  RP11-206L10.8/ENSG00000230092/- 0       -       1       736550  736680  X       0       -       130
#1       736543  741178  RP11-206L10.8/ENSG00000230092/- 0       -       1       736710  736840  X       0       -       130

## Non-Directional
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       134895  135807  E       5
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       135135  135900  E       98
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       135230  136040  X       238
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       136070  136410  X       340
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       136440  136710  X       270
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       136750  137100  X       350
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       137140  137790  X       480
#1       135802  137620  AL627309.1/ENSG00000237683/-    0       -       1       137615  139384  E       5

  my ($intron, $overlapstart, $overlapend, $overlaptype) = $_ =~ /^([^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+\t[^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)\t([^\t]+)/;
  push @overlaps, {'start'=>$overlapstart, 'end'=>$overlapend, 'type'=>$overlaptype};

  if ($lastintron ne $intron) {
    if ($lastintron ne '') {
      processIntron($lastintron, \@overlaps);
      undef @overlaps;
    }
  }
  push @overlaps, {'start'=>$overlapstart, 'end'=>$overlapend, 'type'=>$overlaptype};

  $lastintron = $intron;

}
processIntron($lastintron, \@overlaps);
