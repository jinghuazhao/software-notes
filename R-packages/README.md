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

The currently preferred way to set up is via pkg-config, e.g.,

```bash
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
pkg-config --moderversion jags

R-devel CMD INSTALL --configure-args='--enable-rpath' rjags
```
where `rjags` contains files for the package. Once this is done one can proceed with `install.packages("R2jags")`, etc.
