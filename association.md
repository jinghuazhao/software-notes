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

### VEP

The description is available from http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html.
```bash
git clone https://github.com/Ensembl/ensembl-vep.git

cd ensembl-vep
git pull
git checkout release/92
perl INSTALL.pl
```
The last line requires modules DBI, Build as described in [Overall](../overall).

Lastly, VEP requires .vep directory at $HOME which can be derived from a centrally-installed VEP under Linux,
```bash
ln -s /genetics/ensembl-vep/.vep $HOME/.vep
ln -s /genetics/ensembl-vep/vep $HOME/bin/vep
```
assuming /genetics/ensembl-vep contains the software.

It is slow to get those databases, so one may prefer to get them directly from ftp://ftp.ensembl.org/pub/release-92/variation/VEP/ and unpack into the .vep directory.

We can now execute an example,
```bash
vep -i examples/homo_sapiens_GRCh37.vcf -o out.txt -offline
```

[R-packages](R-packages)

## Pathway analysis

* [DEPICT](DEPICT) (see the [GIANT+Biobank BMI analysis](https://github.com/jinghuazhao/Omics-analysis/tree/master/BMI) ![#f03c15](https://placehold.it/15/f03c15/000000?text=+))
* [PASCAL](PASCAL)

## Medelian randomiszation

* [R-packages](R-packages)

## Polygenic modeling

### [ldsc](ldsc)

### ldetect

It can proceed as indicated
```bash
sudo pip3 install ldetect
```
into /usr/local/lib/python3.6/dist-packages, or
```bash
git clone https://bitbucket.org/nygcresearch/ldetect
cd ldetect
sudo python3 setup.py install
git clone https://bitbucket.org/nygcresearch/ldetect-data
```
A much condensed version of the documentation example is as follows,
```bash
python3 P00_00_partition_chromosome.py example_data/chr2.interpolated_genetic_map.gz 379 example_data/cov_matrix/scripts/chr2_partitions
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr2.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 2:39967768-40067768 | python3 P00_01_calc_covariance.py example_data/chr2.interpolated_genetic_map.gz example_data/eurinds.txt 11418 1e-7 example_data/cov_matrix/chr2/chr2.39967768.40067768.gz
python3 P01_matrix_to_vector_pipeline.py --dataset_path=example_data/cov_matrix/ --name=chr2 --out_fname=example_data/vector/vector-EUR-chr2-39967768-40067768.txt.gz
python3 P02_minima_pipeline.py --input_fname=example_data/vector/vector-EUR-chr2-39967768-40067768.txt.gz  --chr_name=chr2 --dataset_path=example_data/cov_matrix/ --n_snps_bw_bpoints=50 --out_fname=example_data/minima/minima-EUR-chr2-50-39967768-40067768.pickle
python3 P03_extract_bpoints.py --name=chr2 --dataset_path=example_data/cov_matrix/ --subset=fourier_ls --input_pickle_fname=example_data/minima/minima-EUR-chr2-50-39967768-40067768.pickle > example_data/bed/EUR-chr2-50-39967768-40067768.bed
```

### [R-packages](R-packages)

## PheWAS

See wiki resources section of Omics-analysis as well as implementations in [R-packages](R-packages).

## Transcriptome-wide association analysis (TWAS)

* [MetaXcan](MetaXcan) / S-PrediXcan
* [FUSION](FUSION)
* [R-packages](R-packages)

eQTL, single cell sequencing, epigenome-wide association study (EWAS), therapeutic targets
