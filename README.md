viterbi
=======

C++ code for finding shared haplotypes between individuals under the Li and Stephens Model (2003) using the Viterbi algorithm. The program is designed to find long shared haplotypes between a given individual and larger reference panel. 

Usage
-----

To compile, you should just be able to type 'make'. This will create an executable called *viterbi*.

Usage is as follows:

`./viterbi --gzvcf <vcf_1> --gzvcf <vcf_2> --chr <chr> --test-indv <indv_ID> --error <error_rate> > output_file.txt`

+ Parameters
  - --gzvcf : Specify (phased) input VCF file. Can be used more than once, which merged the VCFs internally at biallelic sites with matching alleles. The files must not include overlapping individuals.
  - --chr : Specify the chromosome on which to run the algorithm.
  - --test-indv : Specify the individual for which you want to find the shared haplotypes. (Can be used more than once to specify multiple individuals).
  - --error-rate : Specify the probability of an allele mismatch between shared haplotypes (default: 0.0). 
  - --recomb : Specify the assumed recombination rate in cM/Mb (default: 1.0 cM/Mb).
  - --Ne : Specify the effective population size (default: 10000). 
  - --out : Specify prefix of output files. 

+ Output
The output file has the following columns:
 - HAP : The haplotype in the test individual. The haplotype is suffixed with _1 for the individual's first haplotype, and _2 for the second. 
 - COPY : The haplotype being copied in the reference panel.
 - CHR : Chromosome.
 - START : Haplotype start position.
 - END : Haplotype end position. 
 - N_SNPS : Number of SNPs in the haplotype. 
 - MIN_ALLELE_COUNT : The minimum allele count of alleles present on the haplotype.
 - MIN_ALLELE_POS : The position of the allele with the minimum sample count. 
 - EXACT_HAP_COUNT : The count of the current haplotype in the whole sample.
 - HAP_COUNT_1_PRCT_MISMATCH : The count of the current haplotype in the whole sample, allowing a 1% mismatch. 
  
Notes
-----

This program can use a fair amount of resources. For chr22 on 1000 Genomes Phase 3 data, it needs a few Gb of RAM. For chr1, it needs ~40Gb. A run on chr22 should take about 4 hours on a 3.4Ghz machine, whereas chr1 will take closer to 24 hours.