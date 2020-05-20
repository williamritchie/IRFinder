# IRFinder
Detecting intron retention from RNA-Seq experiments

[User Manual](https://github.com/williamritchie/IRFinder/wiki)

**1.3.0**    
New features:    
1. New `BuildRefFromSTARRef` mode. This allows users to use an existing STAR reference to build IRFinder reference, making the IRFinder result more consistent and comparable with other RNASeq analysis derived from STAR alignment. It also significantly reduces the total preparation time. This new mode also tries to automatically figure out the original FASTA and GTF files used to generate the existing STAR reference. Call `IRFinder -h` for more details.    
2. `BuildRef` and `BuildRefProcess` mode now support `-j` option to parse an integer that changes the default value of `--sjdbOverhang` argument in STAR.    
3. `FASTQ` mode now supports `-y` option to parse a string as extra alignment arguments to STAR when IRFinder quantifies intron retention.    
    
Improvements:    
1. `BAM` mode now outputs a full BAM file in "Unsorted.bam", instead of a BAM file with a trimmed QS column.   
2. IRFinder does not automatically generate "unsorted.frag.bam" to save disk space and to avoid redundancy to "Unsorted.bam". Instead, IRFinder now provides a tool at `bin/TrimBAM4IGV` to generate this kind of trimmed BAM file to facilitate visualization purpose in IGV.     
3. Re-design of standard output information during IRFinder reference preparation. It is easier to recognize occured errors now.    
4. Usage information now can be viewed by `-h` option.     
     
Bug fixes:    
1. The mapability calculation during the IRFinder reference preparation stage has been re-designed. The previous algorithm encountered buffer size issues when dealing with genomes with a huge amount of chromosomes/scaffolds. This has been fixed. Please note, the new algorithm requires `samtools` (>=1.4) executable binary ready in $PATH.    
2. Since Perl 5.28.0, `sort '_mergesort'` is no longer supported. IRFinder now checks the Perl version and uses `sort` functions correspondingly.    
    
**1.2.6**
1. IRFinder now keeps introns with the same effective regions as separate entries in the reference.    
2. IRFinder now automatically checks if the reference preparation stage generates empty reference files, which indicates process failure.    
3. The R object genreated by Differential IR Analysis script now includes an additional slot named "MaxSplice", which represents the maximum splice reads at either end of introns. Each value is the maximum value between Column 17 and 18 in the IR quantification output.    
4. During differential IR analysis, values in "MaxSplice" are now used as the denominators in the GLM, instead of using the values of Column 19 in the IR quantification output. This makes the IR ratio in the differential IR analysis more consistent with the values of Column 20 in the IR quantification output.    
5. User manual has been updated.    
    
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

