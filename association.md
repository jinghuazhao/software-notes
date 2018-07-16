# Association analysis

## --- Single variant analysis ---

### PLINK2

Both [PLINK 1.90 beta](https://www.cog-genomics.org/plink2/) and [PLINK 2.00 alpha](https://www.cog-genomics.org/plink/2.0/) have issue with .grm.bin.N which is shorter than expected for GCTA. The problem is insidious but would prevent chromosome-specific GRMs to be combined.

Nevertheless there is no such problem with its --make-grm-list which allows for the possibility to use --mgrm-list option to combine chromosome-specific GRMs.

Note also the way to use individual's IDs in PLINK2.


### [PyLMM](PyLMM)

## --- HLA imputation ---

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

## --- Finemapping ---

### CAVIAR/eCAVIAR

Installation is made from GitHub in the usual way,
```bash
git clone https://github.com/fhormoz/caviar.git
```

The software requires libgsl which can be installed as follows,
```bash
sudo apt install libgsl-dev
```

Once this is done, one can proceed with the compiling,
```bash
cd caviar
cd CAVIAR-C++
make
```

The references are

**eCAVIAR**

> Hormozdiari F, van de Bunt M, Segrè AV, Li X, Joo JWJ, Bilow M, Sul JH, Sankararaman S, Pasaniuc B, Eskin E (2016). Colocalization of GWAS and eQTL Signals Detects Target Genes. *Am J Hum Genet* 99(6):1245-1260.

**CAVIAR**

> Hormozdiari F, Kostem E, Kang EY, Pasaniuc B, Eskin E (2014). Identifying causal variants at loci with multiple signals of association. *Genetics* 198(2):497-508.

Both are available from https://github.com/fhormoz/caviar.

### [JAM](JAM)

** Setup **

The package is available from https://github.com/pjnewcombe/R2BGLiMS.

Note that JAM requires Java 1.8 so call to Java -jar inside the function needs to
reflect this, not straightforward with `install_github()` from `devtools` but one needs to
clone the package, modify the R source code and then install,
```
git clone https://github.com/pjnewcombe/R2BGLiMS
### change java to java-1.8 in R2BGLiMS/R/R2BGLiMS.R
R CMD INSTALL R2BGLiMS
```

** Compiling **

The information is unavailable from the documentation, but at least can be achieved this with [netbeans](https://netbeans.org/).

## Old notes
```bash
21-7-2017 MRC-Epid JHZ

export BGLiMS=/genetics/bin/BGLiMS
export JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.131-0.b11.el6_9.x86_64

# Jama

cd $BGLiMS
$JAVA_HOME/bin/javac Jama/util/*java
$JAVA_HOME/bin/javac -classpath $BGLiMS Jama/*java
$JAVA_HOME/bin/jar cvf Jama.jar Jama/*class Jama/util/*class

# ssj

export SSJHOME=ssj
export LD_LIBRARY_PATH=$SSJHOME/lib:$LD_LIBRARY_PATH
export CLASSPATH=.:$SSJHOME/lib/ssj.jar:$SSJHOME/lib/colt.jar:$SSJHOME/lib/tcode.jar:$CLASSPATH

# MyClassLibrary

cd $BGLiMS/MyClassLibrary/src
ln -sf $BGLiMS/Jama
$JAVA_HOME/bin/javac Methods/*java
$JAVA_HOME/bin/javac Objects/*java
$JAVA_HOME/bin/jar cvf MyClassLibrary.jar Methods/*class Objects/*class

# BGLiMS.jar

cd $BGLiMS/src
$JAVA_HOME/bin/javac bglims/*java
$JAVA_HOME/bin/jar cvf bglims.jar bglims/*class

cd $BGLiMS
ln -sf $BGLiMS/src/bglims.jar
ln -sf MyClassLibrary/src/MyClassLibrary.jar
# $JAVA_HOME/bin/jar cvf BGLiMS.jar bglims.jar MyClassLibrary.jar
```

## --- Functional annotation ---

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

### [R-packages](R-packages)

## --- Pathway analysis ---

### [DEPICT](DEPICT) (see the [GIANT+Biobank BMI analysis](https://github.com/jinghuazhao/Omics-analysis/tree/master/BMI) ![#f03c15](https://placehold.it/15/f03c15/000000?text=+))

### PASCAL

When there is issue with xianyi-OpenBLAS-v0.2.12-0-g7e4e195.zip shipped with [PASCAL.zip](http://www2.unil.ch/cbg/images/3/3d/PASCAL.zip), 
as described in [vdi.md](https://github.com/jinghuazhao/GDCT/blob/master/vdi.md) or
```{bash}
git clone https://github.com/xianyi/OpenBLAS
```
it is recommended to use the [GitHub version](https://github.com/dlampart/Pascal).

Change to settings.txt is necessary since by default pathway analysis is disabled.

Again we use the BMI summary statistics from GIANT,
```{bash}
wget -qO- http://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz | \
gunzip -c | cut -f1,7 | awk -vFS="\t" -vOFS="\t" '(NR>1)' > BMI.pval

Pascal --pval=BMI.pval

```

## --- Mendelian randomiszation ---

### [R-packages](R-packages)

## --- Polygenic modeling ---

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

### LDpred

It is rather simple to install, e.g.,
```bash
pip install ldpred
```
or
```bash
pip install --user ldpred
```
Nevertheless it is rather laborous to read through the documentation and we try the GitHub version as well.
```bash
git clone https://github.com/bvilhjal/ldpred
cd ldpred
cd ldpred
python test.py
```
rather than a Bash version from scratch which is also possible,
```bash
cd ldpred
export test=ZDNCEO
python coord_genotypes.py --gf=../test_data/LDpred_data_p0.001_train_0 \
                          --vgf=../test_data/LDpred_data_p0.001_test_0 \
                          --ssf=../test_data/LDpred_data_p0.001_ss_0.txt \
                          --N=10000 \
                          --out=$test.coord.hdf5
python LDpred.py --coord=$test.coord.hdf5  \
                 --ld_radius=100 \
                 --local_ld_file_prefix=$test \
                 --PS=0.001 \
                 --N=10000  \
                 --out=$test
python validate.py --vgf=../test_data/LDpred_data_p0.001_test_0 \
                   --rf=$test \
                   --out=$test
```
where
```
 --gf=PLINK_LD_REF_GENOTYPE_FILE,
 --vgf=PLINK_VAL_GENOTYPE_FILE,
 --ssf=SUM_STATS_FILE,
 --NSS_SAMPLE_SIZE, the approximate number of individuals used for calculating the GWAS summary statistics,
 --ld_radiud=LD_RADIUS, the number of SNPs on each side of the focal SNP for which LD should be adjusted,
 --PS=FRACTION_CAUSAL, A list of comma separated (without space) values between 1 and 0, excluding 0. 1 corresponds to the infinitesimal model and will yield results similar to LDpred-inf. Default is --PS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,
 --rf=RESULT_FILE_PREFIX: SNP weights file
```

### [ldsc](ldsc)

### [R-packages](R-packages)

## --- PheWAS ---

See wiki resources section of Omics-analysis as well as implementations in [R-packages](R-packages).

## --- Transcriptome-wide association analysis (TWAS) ---

### [MetaXcan](MetaXcan) / S-PrediXcan

### [FUSION](FUSION)

This section follows http://gusevlab.org/projects/fusion/ and is more compact. To install we do,
```bash
# Software
git clone https://github.com/gusevlab/fusion_twas
cd fusion_twas
# LD reference
wget -qO- https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2 | tar xjvf -
# Gene expression / splicing weights; GTEx weights can be obtained similarly
mkdir WEIGHTS
cd WEIGHTS
for wgt in NTR.BLOOD.RNAARR YFS.BLOOD.RNAARR METSIM.ADIPOSE.RNASEQ CMC.BRAIN.RNASEQ CMC.BRAIN.RNASEQ_SPLICING
do
  wget -qO- https://data.broadinstitute.org/alkesgroup/FUSION/WGT/$wgt.tar.bz2 | tar xfj -
done
# for weight generation only assuming availability of libgfortran.so.3
# wget -qO- https://github.com/genetics-statistics/GEMMA/releases/download/v0.96/gemma.linux.gz | gunzip -c > $HOME/bin/gemma
# ln -s ./ output
```
and add R packages,
```r
library(devtools)
install_github("gabraham/plink2R/plink2R",args="--library=/usr/local/lib/R/site-library/")
install.packages(c('optparse','RColorBrewer'),INSTALL_opts="--library /usr/local/lib/R/site-library/")
# for weight generation
# install.packages('glmnet',INSTALL_opts="--library /usr/local/lib/R/site-library/")
# for joint likelihood mapping
# install_github("cotsapaslab/jlim/jlimR",args="/usr/local/lib/R/site-library/")
# for colocalisation
# install.packages("coloc",INSTALL_opts="/usr/local/lib/R/site-library/")
```
The documentation example for association test, its main use, is then furnished with
```bash
wget https://data.broadinstitute.org/alkesgroup/FUSION/SUM/PGC2.SCZ.sumstats

Rscript FUSION.assoc_test.R \
--sumstats PGC2.SCZ.sumstats \
--weights ./WEIGHTS/NTR.BLOOD.RNAARR.pos \
--weights_dir ./WEIGHTS/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 22 \
--out PGC2.SCZ.22.dat
```

** A useful utility **

```bash
Rscript utils/make_score.R WEIGHTS/CMC.BRAIN.RNASEQ/CMC.MC4R.wgt.RDat > CMC.MC4R.score
plink --bfile genotype-file --score CMC.MC4R.score 1 2 4
```
See additional information from the FUSION documentation.

### [R-packages](R-packages)

eQTL, epigenome-wide association study (EWAS).
