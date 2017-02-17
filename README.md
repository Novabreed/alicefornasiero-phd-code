# phd-code
Collection of code mentioned in Alice Fornasiero PhD thesis.


1.Function clean.geno

Uses thresholds to clean genotypes that were missing in a given number of samples and samples with insuffient genotype information.

Example of input file.
see example for function phase.geno

Parameters:
max.missing.snps: minimum number of polymorphic positions required (default: 0.2)
max.missing.subj: minimum number of genotyped samples required (default: 0.2)



########################
2.Function phase.geno

Define phase of alleles (coded into genotypes) arbitrarily assigned by Stacks

Example of input file.
It contains chromosome, position of RAD locus, joinmap-coded marker type , joinmap-coded genotypes of progenies
CHR POS  MARKER   GENO_SAMPLE1  GENO_SAMPLE2  GENO_SAMPLE3  GENO_SAMPLE4 GENO_SAMPLE5  GENO_SAMPLEN
1   1    <hkxhk>  hk            hh            hh            kk           kk            kk
1   10   <hkxhk>  hk            kk            kk            hh           hh            hh
1   20   <hkxhk>  hk            hh            hh            kk           kk            kk
1   30   <hkxhk>  hk            hh            hh            kk           kk            kk
1   40   <hkxhk>  hk            hh            hh            kk           kk            kk


Parameters:
- threshold.phase: threshold of the cis/trans ratio to assign phase (default: 0.2)
- threshold.homo: minimum number of homozygous genotypes required to assign phase (default: 0.2)
- max.dist.snp: if the information of homozygous genotypes is lacking, it defines the maximum number of previously phased polymorphic positions helping to compute phase of the current position (default: 50)



########################
3.Function define.haplo

Define haplotypes based on homozygous SNPs in the progenies relative to heterozygous SNP positions in the parental variety. 
The function is built to work along with 4 or 5 genotypes.

Example of input file.
It contains chromosome, position, reference of alternate allele and relative genotypes of samples. Genotypes are coded as follows:
0: homozygous for the reference allele
1: heterozygous
2: homozygous for alternate allele

CHR POS  REF ALT GENO_SAMPLE1  GENO_SAMPLE2  GENO_SAMPLE3  GENO_SAMPLE4 (GENO_SAMPLE5)  PARENT_GENO
1   1    A   C   2             0             0             2             (1)            1
1   10   T   G   2             0             0             2             (1)            1
...
1   5000 C  A    0             2             2             0             (0)            1


Parameters:
- samplenames: name of samples analysed. The first name of the list is the parental variety, names of progenies follow.
- remove.het: remove heterozygous genotypes in progenies.
- mystep: define the number of positions within a window to compute the genotype pattern across samples.

