# IRFinder
Detecting intron retention from RNA-Seq experiments

For information on installing and running the software please consult the wiki

**1.2.2:**
1. In GLM-based differential IR comparison, fixed an error caused by duplicated row names when creating DESeq2 object with a version of DESeq2 later than 1.10.

**1.2.1:**
1. Improved the performance of DESeq2-based GLM analysis for differential IR. This new approach should improve the estimation of dispersion. Normal splicing from IRFinder result is now used as a variable in the GLM, instead of using the value of normal splicing as an offset. This approach is adapted from [detection of allele-specific expression](http://rpubs.com/mikelove/ase) from Michael Love. See Wiki page for details.
2. Updated some out-of-date usage information

**1.2.0:**
1. IRFinder is now compatible with GLM-based analysis. This is achieved by passing IRFinder result to DESeq2 using the function in bin/DESeq2Constructor.R. See Wiki page for details  
2. Fixed the conflict with latest version "bedtools complement" that used to cause failure in preparing IRFinder reference  
3. Improved memory usage when passing lines to bedtools genomecov. This is also supposed to benefit reference preparation of those genomes with a lot of chromosomes contigs. Thanks for the smart solution from Andreas @andpet0101.  
4. Specified the gtf file to be downloaded during reference preparation via automatic downloading. Ensembl currently holds several versions of gtf files for the same genome release. This confused IRFinder BuildRefDownload function in the previous version.
5. Added -v option to print out version number.

