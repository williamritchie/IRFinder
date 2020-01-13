#!/bin/bash

export LANG=C
export LC_ALL=C

GTF=$1
CHRLEN=$2
MAPABILITYEXCL=$3
LIBEXEC=$4
OPTNONPOLYA=$5
OPTBLACKLISTEXCL=$6

#GTF=../../REF/Homo_sapiens.GRCh37.75.150/Homo_sapiens.GRCh37.75.limit_chr.gtf
#CHRLEN=../../REF/Homo_sapiens.GRCh37.75.150/STAR/chrNameLength.txt
#MAPABILITYEXCL=MapabilityExclusion.bed.gz
#OPTBLACKLISTEXCL=wgEncodeDacMapabilityConsensusExcludable.bed.gz

NEARGENE="-l 5000 -r 1000 -s"
# Gene has been reversed before being sloped -- effectively adding an exclusion zone 5000 bp downstream of a reverse sense gene (more likely to over-run this direction than to have an early TSS).

echo "Build Ref 1"

"$LIBEXEC/gtf2bed-custom.pl" "$GTF" \
| \
tee \
>(awk 'BEGIN {FS="\t"; OFS="\t"} ($10>1 && $6=="+") {print $1,$2,$3,$4,$5,"-" } ($10>1 && $6=="-") {print $1,$2,$3,$4,$5,"+" }' > tmp.reversed.genes ) \
>(awk 'BEGIN {FS="\t"; OFS="\t"} ($10>1) {print $1,$2,$3,$4,$5,$6 }' > tmp.all.annotations ) \
>(grep -v -e 'intron' | "$LIBEXEC/bed-to-intron+exon.pl" tmp.exons.exclude /dev/null ) \
>(grep -e '/rRNA/' -e '/Mt_rRNA/' | cut -f 1-4,6 | awk 'BEGIN {FS="\t";OFS="\t"} {$4 = "rRNA/" $1 "/" $2 "/" $3 "/" $5 "/" $4; print }' > tmp.ROI.rRNA.bed  ) \
| grep -e '/processed_transcript/' -e '/protein_coding/' | "$LIBEXEC/bed-to-intron+exon.pl" /dev/null tmp.candidate.introns



# grep -e '/rRNA/' -e '/Mt_rRNA/' | cut -f 1-4,6 | sort -S 500M -k1,1 -k2,2n -k3,3n | awk 'BEGIN {FS="\t";OFS="\t"} {$4 = $1 "/" $2 "/" $3 "/" $4; print } '> rRNA-ROI.bed
#>(grep -v -e 'intron' | grep -e '/processed_transcript/' -e '/protein_coding/' | $LIBEXEC/bed-to-intron+exon.pl bed.exons.exclude.anno-coding /dev/null ) \
#>(grep -v -e 'intron' | grep -v -e '/processed_transcript/' -e '/protein_coding/' | $LIBEXEC/bed-to-intron+exon.pl bed.exons.exclude.anno-noncoding /dev/null ) \
#>(awk 'BEGIN {FS="\t"; OFS="\t"} ($10>1 && $4 ~ /protein_coding/) {print $1,$2,$3,$4,$5,$6 }' > tmp.coding.genes ) \

## Stable sort -- lets get the same named introns each time.
sort -t $'\t' -s -S 500M -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 < tmp.candidate.introns | sort -t $'\t' -s -S 500M -k1,1 -k2,2n -k3,3n -k6,6 -u  > introns.unique.bed
#sort -t $'\t' -s -S 500M -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 -u < tmp.candidate.introns > introns.unique.bed

sort -t $'\t' -s -S 500M -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 < tmp.exons.exclude | sort -t $'\t' -s -S 500M -k1,1 -k2,2n -k3,3n -k6,6 -u | bedtools slop -b 5 -i stdin -g "$CHRLEN" | sort -t $'\t' -S 2G -k1,1 -k2,2n -k3,3n -k6,6 > exclude.directional.bed

echo "Build Ref 2"

gzip -dc "$MAPABILITYEXCL" | sort -t $'\t' -S 500M -k1,1 -k2,2n -k3,3n -u | bedtools merge -i stdin -d 9 | awk 'BEGIN {FS="\t"; OFS="\t"} (($3-$2)>=10) { print $1, $2, $3, "M" }' > exclude.omnidirectional.bed

#if [ ! -z "$OPTBLACKLISTEXCL" ]
# if not empty (or null) string.
if [ -n "$OPTBLACKLISTEXCL" ]
then
  if [[ "$OPTBLACKLISTEXCL" == *.gz ]]
  then
    gzip -dc "$OPTBLACKLISTEXCL" | awk 'BEGIN {FS="\t"; OFS="\t"} { if ($1 ~ /^chr/) {sub(/^chr/,"",$1)}; if ($1 == "M") {$1 = "MT"}; print $1,$2,$3,$4 }' >> exclude.omnidirectional.bed
  else
    cat "$OPTBLACKLISTEXCL" | awk 'BEGIN {FS="\t"; OFS="\t"} { if ($1 ~ /^chr/) {sub(/^chr/,"",$1)}; if ($1 == "M") {$1 = "MT"}; print $1,$2,$3,$4 }' >> exclude.omnidirectional.bed
  fi
fi

### BUG in bedtools merge in latest bedtools. Doesn't honour strand (or at least, drops it from the output which is pretty much the same thing)!
echo "Build Ref 3"

function excludeFileDir {
cat \
<( cat exclude.omnidirectional.bed | sort -t $'\t' -S 500M -k1,1 -k2,2n -k3,3n | bedtools merge -i stdin | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, "X", "0", "+"; print $1, $2, $3, "X", "0", "-" }' ) \
<( cat exclude.directional.bed | sort -t $'\t' -S 500M -k1,1 -k2,2n -k3,3n -k6,6 -u | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, "E", "0", $6}' ) \
| \
sort -t $'\t' -S 1G -k1,1 -k2,2n -k3,3n -k6,6
}

function excludeFileNondir {
cat \
<( cat tmp.reversed.genes | sort -t $'\t' -S 500M -k1,1 -k2,2n -k3,3n -k6,6 -u | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, "A", "0", $6}' ) \
<( cat tmp.reversed.genes | bedtools slop $NEARGENE -i stdin -g "$CHRLEN" | sort -t $'\t' -S 2G -k1,1 -k2,2n -k3,3n -k6,6 -u | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, "AE", "0", $6}' ) \
<( cat exclude.omnidirectional.bed | sort -t $'\t' -S 500M -k1,1 -k2,2n -k3,3n -u | bedtools merge -i stdin | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, "X", "0", "+"; print $1, $2, $3, "X", "0", "-" }' ) \
<( cat exclude.directional.bed  | sort -t $'\t' -S 500M -k1,1 -k2,2n -k3,3n -k6,6 -u | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, "E", "0", "+"; print $1, $2, $3, "E", "0", "-" }' ) \
| \
sort -t $'\t' -S 1G -k1,1 -k2,2n -k3,3n -k6,6
}

echo "Build Ref 4"
bedtools intersect -s -sorted -wao -a introns.unique.bed -b <(excludeFileDir) | "$LIBEXEC/IntronExclusion.pl" >(cat > tmp.50) >(cat > tmp.read-continues) | sort -t $'\t' -s -S 500M -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 -u > tmp-dir.IntronCover.bed
echo "Build Ref 5"
bedtools intersect -s -sorted -wao -a introns.unique.bed -b <(excludeFileNondir) | "$LIBEXEC/IntronExclusion.pl" >(cat >> tmp.50) >(cat >> tmp.read-continues) | sort -t $'\t' -s -S 500M -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 -u > tmp-nd.IntronCover.bed
#echo "Build Ref 6"
echo "Build Ref 6"

sort -t $'\t' -S 500M -k1,1 -k2,2n < tmp.read-continues | uniq > ref-read-continues.ref

#awk 'BEGIN {FS = "[\t/]"; OFS = "\t"} {print $1, $8, $9, $6}' < tmp-dir.IntronCover.bed | sort -t $'\t' -S 2G -k1,1 -k2,2n -k3,3n -k4,4 -u > ref-sj.ref
cut -f 1,2,3,6 < introns.unique.bed | sort -t $'\t' -S 2G -k1,1 -k2,2n -k3,3n -k4,4 -u > ref-sj.ref

#echo "Build Ref 8"


# Workaround bug in bedtools merge when dealing with directional data (introduced in latest version of bedtools)
#cat \
#<(awk 'BEGIN {FS="\t"; OFS="\t"} ($6 == "+")' < tmp.coding.genes | sort -t $'\t' -S 2G -k1,1 -k2,2n -k3,3n -k6,6 -u | bedtools merge -i stdin | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, "G+", 0, "+"}' ) \
#<(awk 'BEGIN {FS="\t"; OFS="\t"} ($6 == "-")' < tmp.coding.genes | sort -t $'\t' -S 2G -k1,1 -k2,2n -k3,3n -k6,6 -u | bedtools merge -i stdin | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, "G-", 0, "-"}' ) \
#| sort -t $'\t' -S 2G -k1,1 -k2,2n -k3,3n -k6,6 -u > coding.genes.merged.bed

#echo "Build Ref 9"

#bedtools intersect -S -v -wa -a coding.genes.merged.bed -b <(bedtools slop -b 5000 -i exclude.directional.bed -g "$CHRLEN") | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, "G"$6, $5, "+"; print $1, $2, $3, "G"$6, $5, "-"}' > ref-genes-dir-clean.bed
# awk protein coding -- create description for G+ or G-, but output both directions. The result of these counts indicates whether the sequencing was directional or not.

#cat bed.exons.exclude.anno-coding | awk 'BEGIN {FS="\t"; OFS="\t"} (($3-$2) < 120) {print}' > ref-short-exons.bed 

echo "Build Ref 7"

cat <(awk 'BEGIN {FS="\t"; OFS="\t"} {$4 = "dir/" $4; print}' < tmp-dir.IntronCover.bed) \
  <(awk 'BEGIN {FS="\t"; OFS="\t"} {$4 = "nd/" $4; print}' < tmp-nd.IntronCover.bed) \
  <(awk 'BEGIN {FS="\t"; OFS="\t"} {$4 = "skip"; print}' < tmp.50) \
| sort -t $'\t' -s -S 500M -k1,1 -k2,2n -k3,3n -k6,6 | uniq > ref-cover.bed

echo "Build Ref 8"

#in Bedtools v2.26, chromosome-length file need to be sorted (sort -k1,1) according to chromosome names before being passed to complement function
bedtools slop -b 10000 -g "$CHRLEN" -i tmp.all.annotations \
| sort -t $'\t' -S 500M -k1,1 -k2,2n -k3,3n | bedtools merge -d 1000 -i stdin \
| bedtools complement -i stdin -g <(sort -k1,1 "$CHRLEN") \
| awk 'BEGIN {FS="\t"; OFS="\t"} (length($1)<=2) {$4 = "Intergenic/" $1; print $1, $2, $3, $4}' > intergenic.ROI.bed


echo "Build Ref 9"

if [ -n "$OPTNONPOLYA" ]
then
	if [[ "$OPTNONPOLYA" == *.gz ]]
	then
		bedtools intersect -v -sorted -a <(sort -t $'\t' -S 500M -k1,1 -k2,2n < intergenic.ROI.bed) -b <(cat tmp.ROI.rRNA.bed <(gzip -cd "$OPTNONPOLYA") | cut -f 1-3 | sort -t $'\t' -S 500M -k1,1 -k2,2n) > tmp.ROI.combined.bed
		echo "Build Ref 10a"
		bedtools intersect -v -sorted -a <(sort -t $'\t' -S 500M -k1,1 -k2,2n < tmp.ROI.rRNA.bed) -b <(gzip -cd "$OPTNONPOLYA" | cut -f 1-3 | sort -t $'\t' -S 500M -k1,1 -k2,2n) >> tmp.ROI.combined.bed
		echo "Build Ref 11a"
		gzip -cd < "$OPTNONPOLYA" | awk 'BEGIN {FS="\t"; OFS="\t"} {$4 = "NonPolyA/" $1 "/" $2 "/" $3; print $1, $2, $3, $4}' >> tmp.ROI.combined.bed
		echo "Build Ref 12a"
	else
		bedtools intersect -v -sorted -a <(sort -t $'\t' -S 500M -k1,1 -k2,2n < intergenic.ROI.bed) -b <(cat tmp.ROI.rRNA.bed "$OPTNONPOLYA" | cut -f 1-3 | sort -t $'\t' -S 500M -k1,1 -k2,2n) > tmp.ROI.combined.bed
		echo "Build Ref 10b"
		bedtools intersect -v -sorted -a <(sort -t $'\t' -S 500M -k1,1 -k2,2n < tmp.ROI.rRNA.bed) -b <(cat "$OPTNONPOLYA" | cut -f 1-3 | sort -t $'\t' -S 500M -k1,1 -k2,2n) >> tmp.ROI.combined.bed
		echo "Build Ref 11b"
		awk 'BEGIN {FS="\t"; OFS="\t"} {$4 = "NonPolyA/" $1 "/" $2 "/" $3; print $1, $2, $3, $4}' < "$OPTNONPOLYA" >> tmp.ROI.combined.bed
		echo "Build Ref 12b"
	fi
else
	bedtools intersect -v -sorted -a <(sort -t $'\t' -S 500M -k1,1 -k2,2n < intergenic.ROI.bed) -b <(cat tmp.ROI.rRNA.bed | cut -f 1-3 | sort -t $'\t' -S 500M -k1,1 -k2,2n) > tmp.ROI.combined.bed
	echo "Build Ref 10c"
	cat tmp.ROI.rRNA.bed >> tmp.ROI.combined.bed
	echo "Build Ref 11c"
fi
sort -t $'\t' -S 500M -k1,1 -k2,2n -k3,3n < tmp.ROI.combined.bed > ref-ROI.bed

ENDSTAT=0
if [ ! -s exclude.directional.bed ]
then
    ENDSTAT=1
    echo "Error: exclude.directional.bed is empty."
fi
if [ ! -s exclude.omnidirectional.bed ]
then
    ENDSTAT=1
    echo "Error: exclude.omnidirectional.bed is empty."
fi
if [ ! -s introns.unique.bed ]
then
    ENDSTAT=1
    echo "Error: introns.unique.bed is empty."
fi
if [ ! -s ref-cover.bed ]
then
    ENDSTAT=1
    echo "Error: ref-cover.bed is empty."
fi
if [ ! -s ref-read-continues.ref ]
then
    ENDSTAT=1
    echo "Error: ref-read-continues.ref is empty."
fi
if [ ! -s ref-sj.ref ]
then
    ENDSTAT=1
    echo "Error: ref-sj.ref is empty."
fi

if [ $ENDSTAT = 1]
then
    echo "FAILED"
else
    rm tmp.50 tmp-dir.IntronCover.bed tmp-nd.IntronCover.bed tmp.read-continues tmp.candidate.introns tmp.exons.exclude tmp.all.annotations tmp.reversed.genes tmp.ROI.rRNA.bed tmp.ROI.combined.bed
    echo "COMPLETE"
fi
