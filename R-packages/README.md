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

Again it can be installed with `biocLite("garfield")` and the examination of vignette be done similarly as `coloc`.
