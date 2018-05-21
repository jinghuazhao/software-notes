# R packages

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

### PheWAS

See https://github.com/PheWAS/
