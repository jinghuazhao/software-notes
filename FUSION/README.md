# FUSION

This section is similar to http://gusevlab.org/projects/fusion/ but relatively simple. To install we do,
```bash
# Software
git clone https://github.com/gusevlab/fusion_twas
cd fusion_twas
# LD reference
wget -qO- https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2 | tar xjvf -
# Gene expression / splicing weights; GTEx weights can be obtained similarly
for wgt in NTR.BLOOD.RNAARR YFS.BLOOD.RNAARR METSIM.ADIPOSE.RNASEQ CMC.BRAIN.RNASEQ CMC.BRAIN.RNASEQ_SPLICING
do
  wget -qO- https://data.broadinstitute.org/alkesgroup/FUSION/WGT/$wgt.tar.bz2 | tar xfj -
done
# for weight generation only with complaints about libgfortran.so.3
# wget -qO- https://github.com/genetics-statistics/GEMMA/releases/download/v0.96/gemma.linux.gz | gunzip -c > $HOME/bin/gemma
# ln -s ./ output
```
as well as the R packages,
```r
library(devtools)
install_github("gabraham/plink2R/plink2R",args="--library=/usr/local/lib/R/site-library/")
install.packages(c('optparse','RColorBrewer'),INSTALL_opts="--library /usr/local/lib/R/site-library/")
# for weight generation
# install.packages('glmnet',INSTALL_opts="--library /usr/local/lib/R/site-library/")
```
