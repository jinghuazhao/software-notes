## Issues with DEPICT at GitHub

## Setup

* Information on download and usage is here, https://data.broadinstitute.org/mpg/depict/documentation.html, which has links on [DEPICT_v1_rel194.tar.gz](https://data.broadinstitute.org/mpg/depict/depict_download/bundles/DEPICT_v1_rel194.tar.gz).
It contains 1000Genomes and other data unavailable from [depict_140721.tar.bz2](https://data.broadinstitute.org/mpg/depict/depict_140721.tar.bz2).
* It is preferable to use the source package from GitHub while linking the data directory packaged with DEPICT_v1_rel194.tar.gz above. For instance, 
```{bash}
tar xvfz DEPICT_v1_rel194.tar.gz
export CWD=$(pwd)
git clone https://github.com/perslab/depict
cd depict
mv data data.sav
ln -s $CWD/DEPICT/data
ln -s src/python/depict.py $HOME/bin/depict.py
```
where the package is unpacked into the DEPICT/ directory which contains data/ subdirectory. We then create a symbolic link to data/ from the GitHub version, which does not contain data (to save space on GitHub) while allowing for cutoff_type to be p-values in network analysis for instance.

## PLINK

[PLINK-1.9](https://www.cog-genomics.org/plink2/), with --clump option, has to be used -- [PLINK2](https://www.cog-genomics.org/plink/2.0/) drops the --clump option.

## Python

* call to .sort() should have been .sort_values() in network_plot.py and depict_library.py .

* Is it possible to replicate the Supplementary Figure 9 of the Scott paper? The number of significant pathways seemed to fall short of the FDR<=0.05 criterion. See
 [PW-pipeline](https://github.com/jinghuazhao/PW-pipeline) for how to set up.

## template.cfg

This is from src/python rather than .cfg from example.

## Additional notes

You can examine my [PW-pipeline](https://github.com/jinghuazhao/PW-pipeline) repository on other changes I have made.
