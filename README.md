# IRFinder
Detecting intron retention from RNA-Seq experiments

[User Manual](https://github.com/williamritchie/IRFinder/wiki)

**1.2.5**
1. Headers are now correctly added to output files `IRFinder-IR-dir.txt` and `IRFinder-IR-nondir.txt`.

**1.2.4**
1. In the GLM-based method for differential IR comparison, now the orginal matrix for DESeq2 is now made up by IR depth and correct splicing depth. In the previous versions, the latter one is the sum of splicing depth and IR depth. This change is supposed to give a smoother dispersion estimation across all introns.

**1.2.3:**
1. IRFinder now supports GTF attribution tags `gene_type` and `transcript_type` upon the original requirement for typical Ensembl tags `gene_biotype` and `transcript_biotype`. Either of these two pairs is required to correctly build IRFinder reference.    
    
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

