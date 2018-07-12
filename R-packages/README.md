# R packages

Package installation can be done with `install.packages()` from CRAN and
```r
source("https://bioconductor.org/biocLite.R")
```
from Bioconductor.

## MendelianRandomization

The following are necessary to enable its installation,
```{bash}
sudo apt install curl
sudo apt install libcurl4-openssl-dev
sudo apt install libssl-dev
sudo apt install libgmp-dev
```
and then we have
```{r}
install.packages("MendelianRandomization")
```
The vignette (.R, .Rmd, .pdf) can be seen from /usr/local/lib/R/site-library/MendelianRandomization/doc/Vignette_MR.*

We can now call with
```{bash}
rstudio /usr/local/lib/R/site-library/MendelianRandomization/doc/Vignette_MR.R &
```

# TwoSampleMR

This is standard and furnished as follows,
```r
library(devtools)
install_github('MRCIEU/TwoSampleMR')
```
The following is adapted from 
> Dimou NL, Tsilidis KK (2018). A Primer in Mendelian Randomization Methodology with a Focus on Utilizing Published Summary Association Data in Evangelou E (ed) Genetic Epidemiology-Methods and Protocols. Springer, Chapter 13, pp211-230.
```r
# BMI and lung cancer.

library(TwoSampleMR)
ao <- available_outcomes()
subset(ao,id%in%c(2,966))
exposure_dat <- extract_instruments(ao$id[2])
outcome_dat <- extract_outcome_data(exposure_dat$SNP, 966, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results <- mr(dat)
mr_heterogeneity <- mr_heterogeneity(dat)
mr_pleiotropy_test <- mr_pleiotropy_test(dat)
res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)
p1 <- mr_scatter_plot(mr_results, dat)
p2 <- mr_forest_plot(res_single)
p3 <- mr_leaveoneout_plot(res_loo)
p4 <- mr_funnel_plot(res_single)

library(MendelianRandomization)
MRInputObject <- with(dat, mr_input(bx = beta.exposure, bxse = se.exposure, by = beta.outcome, byse = se.outcome, 
                                    exposure = "Body mass index", outcome = "Lung cancer", snps = SNP))
IVW <- mr_ivw(MRInputObject, model = "default", robust = FALSE, penalized = FALSE, weights = "simple", distribution = "normal", alpha = 0.05)
Egger <- mr_egger(MRInputObject, robust = FALSE, penalized = FALSE, distribution = "normal", alpha = 0.05)
MaxLik <- mr_maxlik(MRInputObject, model = "default", distribution = "normal", alpha = 0.05)
Median <- mr_median(MRInputObject, weighting = "weighted", distribution = "normal", alpha = 0.05, iterations = 10000, seed = 314159265)
MR_all <- mr_allmethods(MRInputObject, method = "all")
p <- mr_plot(MRInputObject, error = TRUE, orientate = FALSE, interactive = TRUE, labels = TRUE, line = "ivw")
```

## BLR

An extensive use is reported in the JSS paper from the [physalia](https://github.com/jinghuazhao/physalia) repository.

## PheWAS

See https://github.com/PheWAS/

## coloc

It requires `snpStats` that can be installed with biocLite().

There is complaint about calling vignette() from Ubuntu; however it is otherwise smooth with help.start().

## garfield

Again it can be installed with `biocLite("garfield")` and vignette be seen similarly to `coloc`.

> GWAS analysis of regulatory or functional information enrichment with LD correction. Briefly, it is a method that leverages GWAS findings with regulatory or 
> functional annotations (primarily from ENCODE and Roadmap epigenomics data) to find features relevant to a phenotype of interest. It performs greedy pruning of 
> GWAS SNPs (LD r2 > 0.1) and then annotates them based on functional information overlap. Next, it quantifies Fold Enrichment (FE) at various GWAS significance 
> cutoffs and assesses them by permutation testing, while matching for minor allele frequency, distance to nearest transcription start site and number of LD 
> proxies (r2 > 0.8).

## rjags

The legacy way to install is 
```bash
R-devel CMD INSTALL --configure-args="--with-jags-prefix=/usr/local --with-jags-libdir=/usr/local/lib --with-jags-includedir=/usr/local/include" rjags
```
but it might not work and the currently preferred way to set up is via pkg-config, e.g.,
```bash
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
pkg-config --moderversion jags

R-devel CMD INSTALL --configure-args='--enable-rpath' rjags
```
where `rjags` contains files for the package. Once this is done one can proceed with `install.packages("R2jags")`, etc.

## sva

The package contains function ```ComBat.R``` from https://www.bu.edu/jlab/wp-assets/ComBat/Download.html as described in the following paper.

> Johnson, WE, Rabinovic, A, and Li, C (2007). Adjusting batch effects in microarray expression data using Empirical Bayes methods. Biostatistics 8(1):118-127

