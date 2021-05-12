# HOSE
HOmogeneity of SNP Effect sizes (HOSE)

This tool can be used to test for the homogeneity of odds ratios (ORs)
between GWAS summary statistics The work is based on the Woolf test for hmogeneity of ORs.
A similar method using the raw data has been used by Liu et al. (2012) Hum Genet.
https://pubmed.ncbi.nlm.nih.gov/21858542/

Three input files are required (for file formats see below):
1) summary statistics for GWAS 1
2) summary statistics for GWAS 2
3) minor allele frequencies in the population

If not provided as part of the summary statistics (1 and 2), then the number of cases and controls for each
GWAS must be provided.

Based on the effect size in the GWAS and the minor allele freqency of the SNP,
the tool computes the allele counts for effect and alternative alleles for cases and controls
in both GWAS. These estimated allel counts will be used in the Woolf test for OR homogeneity to
compute a P-value.

## USAGE:

```
gwas_a_fn <- ""
gwas_b_fn <- ""
maf_fn    <- ""

my_inp  <- harmonize_input(gwas_a_fn, gwas_b_fn, maf_fn)
woolf_p <- HOSE_test(my_inp)

#Alternative (in case sample sizes are missing from sumstat 1):
gwas_a_cases <- 1000
gwas_a_controls <- 10000

gwas_a_sampleinfo <- mk_sample_size_info(gwas_a_cases, gwas_a_controls)
woolf_p <- HOSE_test(my_inp, sample1=gwas_a_sampleinfo)
```

### Adjust p-values for genomic inflation
P-values tend to be inflated. We currently use a simple adjustment based on
the genomic inflation factor.

```
woolf_p_adj <- gc_correct(woolf_p[,1])
```

### Manhattan plot
Quickly create a Manhattan plot.

```
HOSE_man(my_inp, woolf_p_adj)
```

## Implemented tests
There are differnt tests implemented at the moment.
1) woolf: Woolf test (default)
2) BD: Breslow-Day test with correction
3) BD_un: Breslow-Day test w/o correction
3) NMETA: inverse variance meta analysis of the differences of the effect sizes
The methods can be selected in the HOSE_test function using the "method" variable 

## Important
At the moment HOSE expects the files to adhere to some specific formatting 
- All files provide positions according to the same genome build
- Only at most one variant per position (i.e, no triallelic SNPs etc.)
- The effect allele is the same in all files (further scripts to help processing files may be added in the future)


## File formats
File formats (can be zipped with .gz ending)
- GWAS summary statistics (10 columns)
columns contain:
1. Chromosome
2. BP position
3. SNP name
4. Effect allele
5. Reference allele
6. effect size (logOR or 'beta')
7. SE of effect size
8. P-value
9. Sample size cases
10. Sample size controls

```
1       100004726       chr1:100004726  A       G       0.0388  0.0323  0.2293  26421   442271
1       100845815       chr1:100845815  G       A       0.0082  0.0318  0.7959  25252   441303
2       100020742       chr2:100020742  A       G       -0.0001 0.0182  0.9953  33674   449056
```

- GWAS summary statistics (8 columns)
same as above, just missing sample sizes for cases/controls. The number of cases and controls will be provided as aparameter.

- Minor Allele Frequency (e.g., from HRC; 6 columns)
1. Chromosome
2. BP position
3. SNP name
4. Effect allele
5. Reference allele
6. minor allele frequency

```
1       798400  rs10900604      G       A       0.219485
1       2082566 rs2257182       C       T       0.37047
```
