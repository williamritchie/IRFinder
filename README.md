# IRFinder
Detecting intron retention from RNA-Seq experiments

For information on installing and running the software please consult the wiki

**1.2.0:**
1. IRFinder is now compatible with GLM-based analysis. This is achieved by passing IRFinder result to DESeq2 using the function in bin/DESeq2Constructor.R. See Wiki page for details  
2. Fixed the conflict with latest version "bedtools complement” that used to cause failure in preparing IRFinder reference  
3. Improve memory usage when passing lines to bedtools genomecov. This is also supposed to benefit reference preparation of those genomes with a lot of chromosomes contigs. Thanks for the smart solution from Andreas @andpet0101.  
2. Fixed the conflict with latest version "bedtools complement" that used to cause failure in preparing IRFinder reference  
3. Improved memory usage when passing lines to bedtools genomecov. This is also supposed to benefit reference preparation of those genomes with a lot of chromosomes contigs. Thanks for the smart solution from Andreas @andpet0101.  
4. Specified the gtf file to be downloaded during reference preparation via automatic downloading. Ensembl currently holds several versions of gtf files for the same genome release. This confused IRFinder BuildRefDownload function in the previous version.
5. Added -v option to print out version number.

