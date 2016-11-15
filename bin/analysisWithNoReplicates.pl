#!/usr/bin/perl
use strict;
use Data::Dumper;
use List::Util qw(max min);
use FindBin qw($RealBin);
use sort 'stable';


my $winflatExec = "winflat";
if ( -x "$RealBin/util/winflat" ) {
	$winflatExec = "$RealBin/util/winflat";
}else{
	#system('which','winflat', '>/dev/null');
	system('which winflat >/dev/null');
	if ($? != 0) {
		print STDERR "FATAL: winflat utiltiy not found.\nSearched at:\n  $RealBin/util/winflat\n  and on the PATH\n";
		exit 1;
	}
}


sub arrayEqual {
    my ($xref, $yref, $maxCompare) = @_;
    return unless  @$xref == @$yref;

	for (my $i = 0; $i < $maxCompare && $i < scalar @$xref; $i++) {
		return unless $xref->[$i] eq $yref->[$i];
	}
    return 1;
}

sub separatedAB {
    my ($arrayref, $aCount, $bCount) = @_;
	## An array with aCount elements followed by bCount elements.
	## All of the A elements need to be < or > all of the B elements.

	if ($arrayref->[0] == $arrayref->[$aCount]) {
		return 0; #neither > or <.
	}elsif ($arrayref->[0] > $arrayref->[$aCount]) {
		for (my $a = 0; $a < $aCount; $a++) {
			for (my $b = $aCount; $b < $aCount+$bCount; $b++) {
				return 0 if (!($arrayref->[$a] > $arrayref->[$b]));
			}
		}
	}else{
		for (my $a = 0; $a < $aCount; $a++) {
			for (my $b = $aCount; $b < $aCount+$bCount; $b++) {
				return 0 if (!($arrayref->[$a] < $arrayref->[$b]));
			}
		}	
	}
	return 1;
}



my $current = "";
#Filehandles.
my $poolA;
my $poolAname;
my $poolB;
my $poolBname;
my @reps;
my @repsFileNames;
my $repsA = 0;
my $repsB = 0;
my @output;

while (scalar @ARGV) {
  my $param = shift @ARGV;
  if ($param =~ m/^\-/) {
    if ($param eq '-A') {
      $current = 'A';
    }elsif ($param eq '-B') {
      $current = 'B';
    }else{
      print STDERR "Invalid parameter: $param\n";
      exit 1;
    }
  }else{
    if ($current eq "") {
      print STDERR "Invalid parameters. eg: -A pooledA/IR-nondir.txt repA1/IR-nondir.txt repA2/IR-nondir.txt -B pooledB/IR-nondir.txt repB1/IR-nondir.txt repB2/IR-nondir.txt\n";
      exit 1;
    }elsif ($current eq "A") {
      if ($poolA) {
        #open $repsA[scalar @repsA], '<', $param or die "Can't open file $param for reading";
        ## Insert an element into the array @reps, after the last A element.
		splice(@repsFileNames, $repsA, 0, $param);
		splice(@reps, $repsA, 0, undef);
		open $reps[$repsA], '<', $param or die "Can't open file $param for reading";

		$repsA++;
      }else{
        open $poolA, '<', $param or die "Can't open file $param for reading";
        $poolAname = $param;
      }
    }elsif ($current eq "B") {
      if ($poolB) {
        ## Add an element to the very end of the array @reps (ie: after all the A elements, and all the B elements)
		@repsFileNames[scalar @repsFileNames]=$param;
        open $reps[scalar @reps], '<', $param or die "Can't open file $param for reading";
        $repsB++
      }else{
        open $poolB, '<', $param or die "Can't open file $param for reading";
        $poolBname = $param;
      }
    }else{
      print STDERR "error in code\n";
      exit 2;
    }
  }
}

( $poolA ) or die "For condition A, must provide a file.";
( $poolB ) or die "For condition B, must provide a file.";
( $repsA == 0 ) or die "For condition A, must provide a single file only.";
( $repsB == 0 ) or die "For condition B, must provide a single file only.";


#print Dumper(\@repsFileNames);
#print Dumper(\@reps);

print "#Condition A: $poolAname\n";
print "#Condition B: $poolBname\n";


print join("\t",
	"Chr", "Start", "End", "Intron-GeneName/GeneID","-","Direction","ExcludedBases",
	"p-diff","p-increased","p-decreased",
	"A-IRratio","A-IRok","A-IntronCover","A-IntronDepth","A-SplicesMax","A-SplicesExact",
	"B-IRratio","B-IRok","B-IntronCover","B-IntronDepth","B-SplicesMax","B-SplicesExact",
	),"\n";


my $lineNumber = 0;
while(<$poolA>) {
  my $pA = $_;
  chomp $pA;
  my $pB = <$poolB>;
  chomp $pB;
  $lineNumber++;

  my @pA = split /\t/, $pA;
  my @pB = split /\t/, $pB;

  if (!( arrayEqual( \@pA, \@pB, 7) )) {
    print STDERR "FATAL: Files do not list records in the same order with identical number of lines.\n";
    print join("\t", @pA[0 .. 6]), "\n";
    print join("\t", @pB[0 .. 6]), "\n";

    exit 1;
  }

#   ## Loop through replicates, fill an array. (check the ~~ 0..6)
#   my @repsIR;
#   foreach(@reps) {
#   	my @fields = split /\t/, <$_>;
#   	if (!( arrayEqual( \@pA, \@fields, 7) )) {
# 	    print STDERR "FATAL: Files do not list records in the same order with identical number of lines.\n";
# 	    print join("\t", @fields[0 .. 6]), "\n";
# 	    print join("\t", @pA[0 .. 6]), "\n";
# 	    exit 1;
#   	}
#   	push @repsIR, @fields[19];
#   }

  ## Do the maths, are the replicates OK?
#   my $ok = ($pA[20] eq "ok" || $pB[20] eq "ok") && ($pA[8] >= 1 || $pB[8] >= 1) && ($pA[19] >= 0.01 || $pB[19] >= 0.01) && separatedAB(\@repsIR, $repsA, $repsB);

  # No replicates. Still do a form of check -- is this intron interesting?
  my $ok = ($pA[20] eq "-" || $pB[20] eq "-") && ($pA[8] >= 1 || $pB[8] >= 1) && (max($pA[16],$pA[17]) >= 10 || max($pB[16],$pB[17]) >= 10) && ($pA[19] >= 0.01 || $pB[19] >= 0.01);


  my $pValUp = 99;
  my $pValDown = 99;

  if ($ok) {
    ## Check if both are sufficiently expressed (either the intron, or the splices)
    if (( $pA[8] >= 1 || max($pA[16],$pA[17]) >= 10 ) && ( $pB[8] >= 1 || max($pB[16],$pB[17]) >= 10 )) {
	  	## calculate the winflat p-value of the difference (from the pooled IRdepth & junctionDepth).
		#print $lineNumber, "\n";
	    open my $winflat, '-|', $winflatExec, '-xvalue', int($pA[8]+0.5), '-yvalue', int($pB[8]+0.5), '-diff', max($pA[16],$pA[17])+int($pA[8]+0.5), max($pB[16],$pB[17])+int($pB[8]+0.5);
		my @winflat = <$winflat>;
		close $winflat;
	    foreach (@winflat) {chomp; s/^.*\).*= *//; s/\W*$//};
		$pValDown = $winflat[0];
		$pValUp = $winflat[1];
	}else{
		## Properly expressed in only one of the samples. Flag as interesting, but not differential IR.
		$pValUp = 33;
		$pValDown = 33;
	}
  }
  my $pValDiff = min($pValUp, $pValDown);

  if ($ok) {
    # [  ] pushes an array ref onto the array.
    push @output, [@pA[0 .. 6],
    	$pValDiff, $pValUp, $pValDown,
    	$pA[19], $pA[20], $pA[7], $pA[8], max($pA[16],$pA[16]), $pA[18],
    	$pB[19], $pB[20], $pB[7], $pB[8], max($pB[16],$pB[16]), $pB[18]
    	];
  }

  ## Max SJ - 17/18
  ## Exact SJ - 29
  ## IRRatio = 20
  ## ok 21
  ## coverage 8
  ## trimmedMean 9
  ## 
}

foreach ( sort { $a->[7] <=> $b->[7] } @output ) {
	print join("\t", @$_), "\n";
}
#print Dumper (\@output);
