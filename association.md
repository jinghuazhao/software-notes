# Association analysis

## Single variant analysis

### PLINK2

Both [PLINK 1.90 beta](https://www.cog-genomics.org/plink2/) and [PLINK 2.00 alpha](https://www.cog-genomics.org/plink/2.0/) have issue with .grm.bin.N which is shorter than expected for GCTA. The problem is insidious but would prevent chromosome-specific GRMs to be combined.

Nevertheless there is no such problem with its --make-grm-list which allows for the possibility to use --mgrm-list option to combine chromosome-specific GRMs.

Note also the way to use individual's IDs in PLINK2.


### [PyLMM](PyLMM)

## HLA imputation

### HLA*IMP:02

Download source from https://oxfordhla.well.ox.ac.uk/hla/tool/main

```bash
sudo apt install libgd
sudo apt install libgtk-3*
```
then start Perl,
```bash
sudo perl -MCPAN -e shell
```
followed by
```perl
install GD
install Alien::wxWidgets
install Moose
install List::MoreUtils
install Wx::Mini
install Wx::Perl::Packager
```
We can also use cpan, but the installation fails under Ubuntu 18.04.

## Finemapping

* [CAVIAR/eCAVIAR](caviar)
* [JAM](JAM)

## Functional annotation

* [R-packages](R-packages)

## Pathway analysis

* [DEPICT](DEPICT) (see the [GIANT+Biobank BMI analysis](https://github.com/jinghuazhao/Omics-analysis/tree/master/BMI) ![#f03c15](https://placehold.it/15/f03c15/000000?text=+))
* [PASCAL](PASCAL)

## Medelian randomiszation

* [R-packages](R-packages)

## Polygenic modeling

* [ldsc](ldsc)
* [LDpred](ldpred)
* [R-packages](R-packages)

## PheWAS

PheWAS catalog, https://phewascatalog.org/

UK BioBank PheWAS, http://www.ukbiobank.ac.uk/tag/phewas/

See implementations in [R-packages](R-packages).

## Transcriptome-wide association analysis (TWAS)

* [MetaXcan](MetaXcan) / S-PrediXcan
* [FUSION](FUSION)
* [R-packages](R-packages)

eQTL, single cell sequencing, epigenome-wide association study (EWAS), therapeutic targets
