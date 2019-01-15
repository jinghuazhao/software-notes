# Association analysis

## --- Data management ---

### bgenix

As documented, the current version requires gcc 4.7* so we proceed as follows,
```bash
hg clone https://gavinband@bitbucket.org/gavinband/bgen -u master
cd bgen

module load gcc/4.7.2

./waf configure --prefix=/scratch/jhz22
./waf
./waf install

```
See https://bitbucket.org/gavinband/bgen/overview.

## --- Single variant analysis ---

### eigensoft

The PCA software for genomewide data is available from https://www.hsph.harvard.edu/alkes-price/software/ as well as Ubuntu.
```bash
sudo apt install eigensoft
```
The executables are eigenstrat, eigenstratQTL, smarteigenstrat, smartpca, pca, etc.

### GEMMA

To build from source, https://github.com/genetics-statistics/GEMMA, the Makefile needs to change in places with OpenBLAS, /opt/OpenBLAS/.

### PyLMM

The software is rare with its setup for GEI studies accounting for polygenic effects.

**pylmm**

It is necessary to get it going with code in the quick guide,

1. git clone https://github.com/nickFurlotte/pylmm
2. make the following changes,
  * build/scripts-2.7/pylmmGWAS.py, line 207, from `keep = True - v` to `keep = True ^ v`
  * pylmm/lmm.py, line 189, from `if X0 == None:` to `if X0.all == None:`; line 193, from `x = True - np.isnan(Y)` to `x = True ^ np.isnan(Y)`; line 272, from `if X == None: X = self.X0t` to `if X.all == None: X = self.X0t`.
  * build/scripts-2.7/input.py, line 190, from `x = True ^ np.isnan(G)` to `x = True ^ np.isnan(G)`.
3. python setup.py install
4. use the documentation call,
```{bash}
pylmmGWAS.py -v --bfile data/snps.132k.clean.noX --kfile data/snps.132k.clean.noX.pylmm.kin --phenofile data/snps.132k.clean.noX.fake.phenos out.foo
```

**pylmm_zarlab**

Make the following changes to lmm and lmmGWAS similar to pylmm, and then issue `bash run_tests.sh`.

```
EReading SNP input...
Read 1219 individuals from data/snps.132k.clean.noX.fam
Reading kinship...
Read the 1219 x 1219 kinship matrix in 1.139s 
1 number of phenotypes read
Traceback (most recent call last):
  File "scripts/pylmmGWAS.py", line 308, in <module>
    keep = True - v
TypeError: numpy boolean subtract, the `-` operator, is deprecated, use the bitwise_xor, the `^` operator, or the logical_xor function instead.
EReading PLINK input...
Read 1219 individuals from data/snps.132k.clean.noX.fam
Traceback (most recent call last):
  File "scripts/pylmmKinship.py", line 127, in <module>
    K_G = lmm.calculateKinshipIncremental(IN, numSNPs=options.numSNPs,
AttributeError: 'module' object has no attribute 'calculateKinshipIncremental'
EReading PLINK input...
Read 1219 individuals from data/snps.132k.clean.noX.fam
Traceback (most recent call last):
  File "scripts/pylmmKinship.py", line 127, in <module>
    K_G = lmm.calculateKinshipIncremental(IN, numSNPs=options.numSNPs,
AttributeError: 'module' object has no attribute 'calculateKinshipIncremental'

======================================================================
ERROR: test_GWAS (tests.test_lmm.test_lmm)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/tests/test_lmm.py", line 41, in test_GWAS
    TS,PS = lmm.GWAS(Y,snps,K,REML=True,refit=True)
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/pylmm/lmm.py", line 192, in GWAS
    L = LMM(Y, K, Kva, Kve, X0)
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/pylmm/lmm.py", line 301, in __init__
    x = True - np.isnan(Y)
TypeError: numpy boolean subtract, the `-` operator, is deprecated, use the bitwise_xor, the `^` operator, or the logical_xor function instead.

======================================================================
ERROR: test_calculateKinship (tests.test_lmm.test_lmm)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/tests/test_lmm.py", line 25, in test_calculateKinship
    K = lmm.calculateKinship(snps)
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/pylmm/lmm.py", line 135, in calculateKinship
    mn = W[True - np.isnan(W[:, i]), i].mean()
TypeError: numpy boolean subtract, the `-` operator, is deprecated, use the bitwise_xor, the `^` operator, or the logical_xor function instead.

======================================================================
ERROR: test_pylmmGWASScript (tests.test_pylmmGWAS.test_pylmmGWAS)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/tests/test_pylmmGWAS.py", line 24, in test_pylmmGWASScript
    with (open(self._outputFile, 'r')) as ansFile:
IOError: [Errno 2] No such file or directory: 'data/pylmmGWASTestOutput'

======================================================================
ERROR: test_pylmmKinshipScript1 (tests.test_pylmmKinship.test_pylmmKinship)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/tests/test_pylmmKinship.py", line 24, in test_pylmmKinshipScript1
    K = np.fromfile(open(self._outputFile, 'r'), sep=" ")
IOError: [Errno 2] No such file or directory: 'data/pylmmKinshipTestOutput'

======================================================================
ERROR: test_pylmmKinshipScript2 (tests.test_pylmmKinship.test_pylmmKinship)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/tests/test_pylmmKinship.py", line 38, in test_pylmmKinshipScript2
    K = np.fromfile(open(self._outputFile, 'r'), sep=" ")
IOError: [Errno 2] No such file or directory: 'data/pylmmKinshipTestOutput'

----------------------------------------------------------------------
Ran 5 tests in 2.776s

FAILED (errors=5)
```
We can have a test of GxE analysis as this,
```{bash}
sudo python setup.py install
cd pylmm
python pylmm_GXE.py
```
In general, we can see options for GxE analysis from command `pylmmGWAS.py` under bash.

### METAL

Note that at least cmake 3.1 is required for the latest from GitHub, https://github.com/statgen/METAL, 

```bash
wget -qO- https://github.com/statgen/METAL/archive/2018-08-28.tar.gz | \
tar xvfz -
cd METAL-2018-08-28
mkdir build && cd build
cmake ..
make
make test
make install
```
One can use `ccmake .` to change the prefix for installation but this does not appear to work and the executable is bin/metal.

As with distribution 2011-03-25, http://csg.sph.umich.edu/abecasis/Metal/download/, options CUSTOMVARIABLE uses an output format of %g, leading to scientific notation of position, which
is undesirable and we modify metal/Main.cpp from
```c
for  (int j = 0; j < customVariables.Length(); j++)
            fprintf(f, "\t%g", custom[j][marker]);

```
to
```c
for  (int j = 0; j < customVariables.Length(); j++)
            fprintf(f, "\t%-.15g", custom[j][marker]);

```
which is left-aligned with 15 places with %g though largely 11 is enough. The change can 
be tested by adding the following lines to examples/GlucoseExample/meta.txt.
```
CUSTOMVARIABLE CHR
LABEL CHR as CHR
CUSTOMVARIABLE POS
LABEL POS as POS
CUSTOMVARIABLE N
LABEL N as N

CHROMOSOMELABEL CHR
POSITIONLABEL POS
TRACKPOSITIONS ON
```
Nevertheless the CHR and POS thus retained are sums of individual studies involved for particular positions so their real values
can be recovered from these divided by the number of - and + from the Direction column. In the case of N, the sum is just what
we want. To wrap up, our testing code is as follows,
```bash
### illustration as in examples/GlucoseExample/metal.tbl of TRACKPOSITIONS and change in source
### the results are also sorted in accordance with METAL documentation

cd examples/GlucoseExample
(
echo CUSTOMVARIABLE CHR
echo LABEL CHR as CHR
echo CUSTOMVARIABLE POS
echo LABEL POS as POS
echo CUSTOMVARIABLE N
echo LABEL N as N

echo CHROMOSOMELABEL CHR
echo POSITIONLABEL POS
echo TRACKPOSITIONS ON

echo OUTFILE metal- .tbl

awk '!/\#/' metal.txt 
) > metal.metal

metal metal.metal
(
head -1 metal-1.tbl
awk 'NR>1' metal-1.tbl | \
awk '
{
   FS=OFS="\t"
   direction=$9
   gsub(/\?/,"",direction)
   n=length(direction)
   $10=$10/n
   $11=$11/n
};1' | \
sort -k10,10n -k11,11n
) > metal.tbl
rm metal-1.tbl
mv metal-1.tbl.info metal.tbl.info

cd -
```

---

### mtag

https://github.com/omeed-maghzian/mtag

### seqMeta

seqMeta: Meta-Analysis of Region-Based Tests of Rare DNA Variants, https://cran.r-project.org/web/packages/seqMeta/index.html.

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

---

## --- Finemapping ---

### CAVIAR/eCAVIAR

Installation is made from GitHub in the usual way,
```bash
git clone https://github.com/fhormoz/caviar.git
```

The software requires libgsl and liblapack which can be installed as follows,
```bash
sudo apt install liblapack-dev
sudo apt install libgsl-dev
```

Once this is done, one can proceed with the compiling,
```bash
cd caviar
cd CAVIAR-C++
make
```
It may be necessary to alter Makefile to point to appropriate -I -L for lapack, for instance.

The references are

**eCAVIAR**

> Hormozdiari F, van de Bunt M, Segrè AV, Li X, Joo JWJ, Bilow M, Sul JH, Sankararaman S, Pasaniuc B, Eskin E (2016). Colocalization of GWAS and eQTL Signals Detects Target Genes. *Am J Hum Genet* 99(6):1245-1260.

**CAVIAR**

> Hormozdiari F, Kostem E, Kang EY, Pasaniuc B, Eskin E (2014). Identifying causal variants at loci with multiple signals of association. *Genetics* 198(2):497-508.

Both are available from https://github.com/fhormoz/caviar.

### CAVIARBF

```bash
wget https://bitbucket.org/Wenan/caviarbf/get/7e428645be5e.zip
unzip 7e428645be5e.zip
cd Wenan-caviarbf-7e428645be5e
make
ln -sf $PWD/caviarbf $HOME/bin/caviarbf
ln -sf $PWD/model_search $HOME/bin/model_search
./install_r_package.sh 
cd caviarbf-r-package
R --no-save <<END
install.packages("glmnet")
END
R CMD INSTALL caviarbf_0.2.1.tar.gz 
./test.sh
cd -
./test.sh
```

### JAM

**Setup**

The package is available from https://github.com/pjnewcombe/R2BGLiMS.

Note that JAM requires Java 1.8 so call to Java -jar inside the function needs to
reflect this, not straightforward with `install_github()` from `devtools` but one needs to
clone the package, modify the R source code and then install,
```
git clone https://github.com/pjnewcombe/R2BGLiMS
### in case you have java-1.6 you will need to change to java-1.8 in R2BGLiMS/R/R2BGLiMS.R
### and possibly add other options, e.g.,
### /usr/lib/jvm/java-8-oracle/bin/java -Xmx4G
### sed -i 's|\"java|\"/usr/lib/jvm/java-8-oracle/bin/java -Xmx4G||g' R2BGLiMS/R/R2BGLiMS.R
### sudo R CMD INSTALL R2BGLiMS -l $R_LIBS
### if R_LIBS is not set, the default can be used, e.g., $HOME/R
R CMD INSTALL R2BGLiMS
```
As shown, it might well be necessary to add options to the Java command-line.

**Compiling**

The information is unavailable from the documentation, but can be achieved this with [netbeans](https://netbeans.org/) or in steps.
```bash
# 21-7-2017 MRC-Epid JHZ

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

---

## --- Functional annotation ---

### fgwas

```bash
git clone https://github.com/joepickrell/fgwas
cd fgwas
# sed -i 's/1.14/1.15/g' configure
./configure --prefix=/scratch/jhz22
make
sudo make install
src/fgwas -i test_data/test_LDL.fgwas_in.gz -w ens_coding_exon
git clone https://github.com/joepickrell/1000-genomes
git clone https://github.com/joepickrell/1000-genomes-genetic-maps
```
In case boots is unavailable or not up-to-date, it is necessary to install, e.g.,
```bash
wget https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz
tar xvfz boost_1_69_0.tar.gz
cd boost_1_60_0
./bootstrap.sh --prefix=/scratch/jhz22 --exec-prefix=/scratch/jhz22
./b2 install
```
Nevertheless it is also required to have other dependencies in place.

### VEP

The description is available from http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html.
```bash
git clone https://github.com/Ensembl/ensembl-vep.git

cd ensembl-vep
git pull
git checkout release/92
perl INSTALL.pl
```
The last line requires modules DBI, Build as described in the [LANGUAGES](https://github.com/jinghuazhao/Computational-Statistics/blob/master/LANGUAGES.md) section of [Computational-Statistics](https://github.com/jinghuazhao/Computational-Statistics).

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

### R-packages

See R-packages section.

---

## --- Pathway analysis ---

### DEPICT

See the [GIANT+Biobank BMI analysis](https://github.com/jinghuazhao/Omics-analysis/tree/master/BMI) ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)

**Installation and documentation example**

* The official site, https://data.broadinstitute.org/mpg/depict/documentation.html has links on [DEPICT_v1_rel194.tar.gz](https://data.broadinstitute.org/mpg/depict/depict_download/bundles/DEPICT_v1_rel194.tar.gz),
which contains 1000Genomes and other data unavailable from [depict_140721.tar.bz2](https://data.broadinstitute.org/mpg/depict/depict_140721.tar.bz2).
```bash
wget https://data.broadinstitute.org/mpg/depict/depict_download/bundles/DEPICT_v1_rel194.tar.gz
tar xvfz DEPICT_v1_rel194.tar.gz
export CWD=$(pwd)
```
where the package is unpacked into the DEPICT/ directory containing the data/ subdirectory. We also note down current working directory with `CWD`.

* The source package from GitHub has more features such as cutoff_type to be p-values in network analysis; the code
```{bash}
git clone https://github.com/perslab/depict
cd depict
wget https://data.broadinstitute.org/mpg/depict/depict_download/collections/ld0.5_collection_1000genomespilot_depict_150429.txt.gz
mkdir -p data/collections
mv ld0.5* data/collections
sed 's|/cvar/jhlab/tp/DEPICT|/home/jhz22/Downloads/depict|g;s|label_for_output_files: ldl_teslovich_nature2010|label_for_output_files: test|g; s|/cvar/jhlab/tp/tools/plink/plink-1.09-Sep2015-x86_64/plink|/home/jhz22/bin/plink|g' example/ldl_teslovich_nature2010.cfg > test.cfg
src/python/depict.py test.cfg
```
adds `ld0.5_collection_1000genomespilot_depict_150429.txt.gz` and produces results prefixed with `test_` using the LDL data.

* Since the documentation example above does not give the full results, data directory packaged with DEPICT_v1_rel194.tar.gz above is called to remedy with a minor change,
```bash
mv data data.sav
ln -s $CWD/DEPICT/data
```
to `test.cfg` for a re-run.
```bash
sed -i 's|data/reconstituted_genesets/reconstituted_genesets_example.txt|data/reconstituted_genesets/reconstituted_genesets_150901.binary|g' test.cfg
src/python/depict.py test.cfg
```

* PLINK. [PLINK-1.9](https://www.cog-genomics.org/plink2/), with --clump option, has to be used rather than [PLINK2](https://www.cog-genomics.org/plink/2.0/) since itdrops the --clump option.

* NB template.cfg is from src/python rather than .cfg from example.

* Python 2.7.*. After installation, the following change is needed: from .sort() to .sort_values() in network_plot.py and depict_library.py. It is necessary to download [additional files](https://data.broadinstitute.org/mpg/depict/depict_download/) for network analysis -- in my case, downloads via Firefox do not work and I used `wget` instead.

* To explore possibility to replicate the Supplementary Figure 9 of the Scott paper -- the number of significant pathways seemed to fall short of the FDR<=0.05 criterion, see
[SUMSTATS](https://github.com/jinghuazhao/SUMSTATS) for setup.

* Under Windows, `gzip.exe` is also required at the working directory or %path% plus some changes over directory specifications. We can then execute
```
python depict.py BMI.cfg
```

* For tissue plot, one can use pdftopng from XpdfReader (or convert/magick from ImageMagick) to obtain .png files to be incorporated into Excel workbook. For network plot, the python package scikit-learn is required.
```bash
sudo pip install scikit-learn
```

**Recompile**

This may be necessary for large collection of significant variants, e.g., GIANT+UKB height summary statistics (height_loci.txt has 2,184 lines including header).

Start netbeans and open project from depict/src/java, fixing links to colt.jar, commons-math-2.0.jar, Jama-1.0.2.jar, jsci-core.jar, JSciCore.jar, jsc.jar from lib/.

**Additional notes**

[PW-pipeline](https://github.com/jinghuazhao/PW-pipeline) puts together many changes and is streamlined with other software.


### MAGMA

A generic setup is available from [PW-pipeline](https://github.com/jinghuazhao/PW-pipeline), while section CAD of the
[Omics-analysis](https://github.com/jinghuazhao/Omics-analysis) repository provides a much simplified version.

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

### VEGAS2

It is relatively slow with web interface https://vegas2.qimrberghofer.edu.au, so we would like to try the command-line counterpart. Make sure R packages `corpcor` and `mvtnorm` are available, then proceed with
```bash
# driver download
wget https://vegas2.qimrberghofer.edu.au/vegas2v2
# documentation example -- unzip does not accept input from console so we do in two steps
wget https://vegas2.qimrberghofer.edu.au/VEGAS2v2example.zip
unzip -j VEGAS2v2example.zip
# gene-based association
perl vegas2v2 -G -snpandp example.txt -custom $PWD/example -glist example.glist -genelist example.genelist -out example
# pathway-based association
awk '(NR>1){OFS="\t";gsub(/"/,"",$0);print $2,$8}' example.out > Example.geneandp
vegas2v2 -P -geneandp Example.geneandp -glist example.glist -geneandpath Example.vegas2pathSYM -out Example
# further setup
wget https://vegas2.qimrberghofer.edu.au/biosystems20160324.vegas2pathSYM
wget https://vegas2.qimrberghofer.edu.au/glist-hg19
wget -qO- https://vegas2.qimrberghofer.edu.au/g1000p3_EUR.tar.gz | tar xvfz -
```
Somehow the binary files following `-custom` option needs to be absolute path.

The last line downloads and unpacks the LD reference data for European (EUR) population; other options include AFR, AMR, EAS, SAS.

---

## --- Mendelian randomiszation ---

See R-packages section.

## --- Polygenic modeling ---

### HESS

HESS (Heritability Estimation from Summary Statistics) is now available from https://github.com/huwenboshi/hess and has a web page at

https://huwenboshi.github.io/hess-0.5/#hess

Some popular Python packages such as pandas as well as PySnpTools are required, e.g., `pip install pandas` or `python -m pip install pysnptools` and it is doable for the latter with source at https://github.com/MicrosoftGenomics/PySnpTools.

For the GIANT height data, we had success with the following script,
```bash
#!/bin/bash

export HEIGHT=https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz

wget -qO- $HEIGHT | \
awk 'NR>1' | \
sort -k1,1 | \
join -13 -21 snp150.txt - | \
awk '($9!="X" && $9!="Y" && $9!="Un"){if(NR==1) print "SNP CHR BP A1 A2 Z N"; else print $1,$2,$3,$4,$5,$7/$8,$10}' > height.tsv.gz

#  SNP - rs ID of the SNP (e.g. rs62442).
#  CHR - Chromosome number of the SNP. This should be a number between 1 and 22.
#  BP - Base pair position of the SNP.
#  A1 - Effect allele of the SNP. The sign of the Z-score is with respect to this allele.
#  A2 - The other allele of the SNP.
#  Z - The Z-score of the SNP.
#  N - Sample size of the SNP.

for chrom in $(seq 22)
do
    python hess.py \
        --local-hsqg height \
        --chrom $chrom \
        --bfile 1kg_eur_1pct/1kg_eur_1pct_chr${chrom} \
        --partition nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-chr${chrom}.bed \
        --out step1
done
python hess.py --prefix step1 --out step2
```
where snp150.txt from UCSC is described at the SUMSTATS repository, https://github.com/jinghuazhao/SUMSTATS.

### ldetect

It can proceed as indicated
```bash
sudo pip3 install ldetect
```
into /usr/local/lib/python3.6/dist-packages, or
```bash
git clone https://bitbucket.org/nygcresearch/ldetect
cd ldetect
# for super user
# sudo python3 setup.py install
python3 setup.py install --user
git clone https://bitbucket.org/nygcresearch/ldetect-data
cd ldetect-data
for pop in AFR ASN EUR
do awk '
{
  OFS="\t"
  if (NR==1) print "#chrom", "Start", "End", "Region"
  else  print $1, $2, $3, "region" NR-1
}' $pop/fourier_ls-all.bed > $pop.bed
done
cd -
```
Population-specific approximately independent LD blocks are given. A much condensed version of the documentation example is as follows,
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

### ldsc

**Partitioned heritability**

The [wiki documentation](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability) script really should be as follows,
```bash

export GIANT_BMI=GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt

setup()
{
  if [ ! -f $GIANT_BMI ]; then
     wget http://portals.broadinstitute.org/collaboration/giant/images/b/b7/$GIANT_BMI.gz
     gunzip $GIANT_BMI.gz
  fi
  if [ ! -f w_hm3.snplist ]; then
     wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
     bzip2 -d w_hm3.snplist.bz2
  fi
  python munge_sumstats.py --sumstats $GIANT_BMI --merge-alleles w_hm3.snplist --out BMI --a1-inc
}

# It fails to munge, so we use brute-force. For P-to-Z implementation in C/C++, see
# https://stackoverflow.com/questions/27830995/inverse-cumulative-distribution-function-in-c
# https://stackoverflow.com/questions/22834998/what-reference-should-i-use-to-use-erf-erfc-function
# It turned out that an older version of pandas is required here, see the GIANT+UKB BMI example

awk 'NR>1' $GIANT_BMI > 1
awk 'NR>1' w_hm3.snplist | sort -k1,1 | join -j1 1 - | awk -f CLEAN_ZSCORES.awk > BMI.sumstats
R --vanilla -q <<END
BMI <- read.table("BMI.sumstats",col.names=c("SNP","A1","A2","Z","N"))
BMI <- within(BMI, {Z=sign(Z)*qnorm(abs(Z)/2)})
z <- gzfile("BMI.sumstats.gz","w")
write.table(BMI,file=z,quote=FALSE,row.names=FALSE)
close(z)
END
```
where we use [CLEAN_ZSCORES.awk](files/CLEAN_ZSCORES.awk) to align SNPs between sumstats and reference.

Now the partition heritability and cell-type group analysis proceed as follows,
```bash
python ldsc.py --h2 BMI.sumstats.gz\
        --ref-ld-chr baseline_v1.1/baseline.\
        --w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.\
        --overlap-annot\
        --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC.\
        --out BMI_baseline

python ldsc.py --h2 BMI.sumstats.gz\
        --w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.\
        --ref-ld-chr 1000G_Phase3_cell_type_groups/cell_type_group.3.,baseline_v1.1/baseline.\
        --overlap-annot\
        --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC.\
        --out BMI_CNS\
        --print-coefficients
```
NB it is assumed that [all the required data](https://data.broadinstitute.org/alkesgroup/LDSCORE/) have been made available.

**Test**

The mysterious commands shown in the wiki documentation are actually realised after this,
```
sudo apt install python-nose
```

### PLINK2

Both [PLINK 1.90 beta](https://www.cog-genomics.org/plink2/) and [PLINK 2.00 alpha](https://www.cog-genomics.org/plink/2.0/) have issue with .grm.bin.N which is shorter than expected for GCTA. The problem is insidious but would prevent chromosome-specific GRMs to be combined.

Nevertheless there is no such problem with its --make-grm-list which allows for the possibility to use --mgrm-list option to combine chromosome-specific GRMs.

Note also the way to use individual's IDs in PLINK2.

### R-packages

See R-packages section.

---

## --- PheWAS ---

See wiki resources section of Omics-analysis as well as implementations in R-packages section.

---

## --- Transcriptome-wide association analysis (TWAS) ---

### MetaXcan / S-PrediXcan

**Issues with more recent version of Python 2.7**

The issue was raised to the [MetaXcan GitHub repository](https://github.com/hakyimlab/MetaXcan) for

* Ubuntu 18.04 and Python 2.7.15r1
* Fedora 27 and Python 2.7.15

such that logging.getLogger() was not found.

This was due to confusion between Logging and Python module logging, and fixed by renaming Logging.py to myLogging.py and then adjusting the call from Logging to 
myLogging in M03_betas.py, M04_zscores.py and MetaXcan.py, etc.

It looks the recent version of Python is stricter, which is somewhat expected as with most other compilers. Similar issues were raised while maintaining R packages
for complaints from g++ 8.xx (to be shipped with Fedora 28) which is otherwise OK with g++ 7.x.x.

**Use of the latest databases**

While it is possible to use the web interface, https://cloud.hakyimlab.org/user_main, to achieve greater flexibility, the latest databases can be downloaded locally
from [PredictDB Data Repository](http://predictdb.org/).

For instance with [GTEx-V7_HapMap-2017-11-29.tar.gz](https://s3.amazonaws.com/predictdb2/GTEx-V7_HapMap-2017-11-29.tar.gz), we can do the following steps,
```{bash}
mkdir GTEx-V7_HapMap-2017-11-29
cd GTEx-V7_HapMap-2017-11-29
wget https://s3.amazonaws.com/predictdb2/GTEx-V7_HapMap-2017-11-29.tar.gz
tar xvfz GTEx-V7_HapMap-2017-11-29.tar.gz
```
and adjust for the documentation example 
```{python}
./MetaXcan.py \
--model_db_path data/DGN-WB_0.5.db \
--covariance data/covariance.DGN-WB_0.5.txt.gz \
--gwas_folder data/GWAS \
--gwas_file_pattern ".*gz" \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column BETA \
--pvalue_column P \
--output_file results/test.csv
```
as follows,
```{python}
./MetaXcan.py \
--model_db_path /home/jhz22/D/genetics/hakyimlab/ftp/GTEx-V7_HapMap-2017-11-29/gtex_v7_Brain_Amygdala_imputed_europeans_tw_0.5_signif.db \
--covariance /home/jhz22/D/genetics/hakyimlab/ftp/GTEx-V7_HapMap-2017-11-29/gtex_v7_Brain_Amygdala_imputed_eur_covariances.txt.gz \
--gwas_folder data/GWAS \
--gwas_file_pattern ".*gz" \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column BETA \
--pvalue_column P \
--output_file results/V7.csv
```

**Examining weights and related information**

[PredictDB FAQs](http://predictdb.org/FAQ.html) point to [a utility in PrediXcan](https://github.com/hakyimlab/PrediXcan/blob/master/Software/query-db.Rmd) for
query, however it is handy to use sqlite3 directory as has been demonstrated in my [TWAS-pipeline](https://github.com/jinghuazhao/TWAS-pipeline). In this case,
we can create a utility, called `query-db.sql` here, 
```{sqlite3}
.tables
.separator "\t"
.header on
.output weights.txt
select * from weights;
.output extra.txt
select * from extra;
```
used as follows,
```{bash}
sqlite3 gtex_v7_Brain_Amygdala_imputed_europeans_tw_0.5_signif.db < query-db.sql
```
and the weights and extra information are available from files weights.txt and extra.txt, respectively.

### FUSION

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

**A useful utility**

```bash
Rscript utils/make_score.R WEIGHTS/CMC.BRAIN.RNASEQ/CMC.MC4R.wgt.RDat > CMC.MC4R.score
plink --bfile genotype-file --score CMC.MC4R.score 1 2 4
```
See additional information from the FUSION documentation.

### ExPecto

Software for predicting expression effects of human genome variants ab initio from sequence.
```bash
git clone https://github.com/FunctionLab/ExPecto
cd ExPecto
sudo pip install -r requirements.txt
sh download_resources.h
tar fxz resources.tar.gz
python chromatin.py example/example.vcf
python predict.py --coorFile example/example.vcf --geneFile example/example.vcf.bed.sorted.bed.closestgene --snpEffectFilePattern example/example.vcf.shift_SHIFT.diff.h5 --modelList resources/modellist --output output.csv
python train.py --expFile resources/geneanno.exp.csv --targetIndex 1 --output model.adipose
```

### INFERNO

Short for (INFERring the molecular mechanisms of NOncoding genetic variants, it is available from https://bitbucket.org/wanglab-upenn/INFERNO 
and it also has a web interface, http://inferno.lisanwanglab.org/index.php.

### R packages

pSI, available from CRAN and http://genetics.wustl.edu/jdlab/psi_package/ with supplementary data http://genetics.wustl.edu/jdlab/files/2014/01/pSI.data_1.0.tar_.gz.

---

## R-packages

See https://github.com/jinghuazhao/Computational-Statistics for general information.

**Bioconductor**

This includes Biobase, BSGenome, edgeR, limma, Rsubread, STRINGdb.

**CRAN**

This includes DCGL.

**GSMR**

It refers to Generalised Summary-data-based Mendelian Randomisation, http://cnsgenomics.com/software/gsmr/, available both as part of GCTa and R package.
```r
install.packages("http://cnsgenomics.com/software/gsmr/static/gsmr_1.0.6.tar.gz",repos=NULL,type="source")
```
with test data, http://cnsgenomics.com/software/gsmr/static/test_data.zip.

**MendelianRandomization**

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
We can use the 97 SNPs from GIANT as described in [SUMSTATS](https://github.com/jinghuazhao/SUMSTATS) as two subsets and obtain association 
informaiton in _PhenoScanner_GWAS.csv using [PhenoScanner](http://www.phenoscanner.medschl.cam.ac.uk/) as well as the `extract.pheno.csv()`
to build a MR analysis for BMI-T2D, say.


**TwoSampleMR**

This is standard and furnished as follows,
```r
library(devtools)
install_github('MRCIEU/TwoSampleMR')
```
The following is adapted from 
> Dimou NL, Tsilidis KK (2018). A Primer in Mendelian Randomization Methodology with a Focus on Utilizing Published Summary Association Data in Evangelou E (ed) Genetic Epidemiology-Methods and Protocols. Springer, Chapter 13, pp211-230.
```r
# BMI and T2D.

library(TwoSampleMR)
ao <- available_outcomes()
subset(ao,id%in%c(2,24))
ao2 <- subset(ao,id==2)
exposure_dat <- extract_instruments(ao2$id)
outcome_dat <- extract_outcome_data(exposure_dat$SNP, 24, proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
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
                                    exposure = "Body mass index", outcome = "Type 2 diabetes", snps = SNP))
IVW <- mr_ivw(MRInputObject, model = "default", robust = FALSE, penalized = FALSE, weights = "simple", distribution = "normal", alpha = 0.05)
Egger <- mr_egger(MRInputObject, robust = FALSE, penalized = FALSE, distribution = "normal", alpha = 0.05)
MaxLik <- mr_maxlik(MRInputObject, model = "default", distribution = "normal", alpha = 0.05)
Median <- mr_median(MRInputObject, weighting = "weighted", distribution = "normal", alpha = 0.05, iterations = 10000, seed = 314159265)
MR_all <- mr_allmethods(MRInputObject, method = "all")
p <- mr_plot(MRInputObject, error = TRUE, orientate = FALSE, interactive = TRUE, labels = TRUE, line = "ivw")
pdf("BMI-T2D.pdf")
p1
p2
p3
p4
p
dev.off()
```
[ACE-APOE-CRP.R](files/ACE-APOE-CRP.R) illustrates the use of MRInstruments, linking some established proteins.

**BLR**

An extensive use is reported in the JSS paper from the [Mixed-Models](https://github.com/jinghuazhao/Mixed-Models) repository.

**PheWAS**

See https://github.com/PheWAS/

**coloc**

It requires `snpStats` that can be installed with biocLite().

There is complaint about calling vignette() from Ubuntu; however it is otherwise smooth with `help.start()`.

Here we run examples modified from the documentation,
```r
# coloc, large (>0.05) p.value.chisquare indicates traits are compatible with colocalisation
# https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html
set.seed(1)
X1 <- matrix(rbinom(1200,1,0.4),ncol=2)
X2 <- matrix(rbinom(1000,1,0.6),ncol=2)
colnames(X1) <- colnames(X2) <- c("f1","f2")
Y1 <- rnorm(600,apply(X1,1,sum),2)
Y2 <- rnorm(500,2*apply(X2,1,sum),5)
summary(lm1 <- lm(Y1~f1+f2,data=as.data.frame(X1)))
summary(lm2 <- lm(Y2~f1+f2,data=as.data.frame(X2)))
require(coloc)
## intuitive test for proportionality
ct <- coloc.test(lm1,lm2, plots.extra=list(x=c("eta","theta"), y=c("lhood","lhood")))
summary(ct)
b1 <- coef(lm1)
b2 <- coef(lm2)
v1 <- vcov(lm1)
v2 <- vcov(lm2)
coloc.test.summary(b1,b2,v1,v2)
# some Bayesian flavour
ct.bayes <- coloc.test(lm1,lm2, plots.extra=list(x=c("eta","theta"), y=c("lhood","lhood")),bayes=TRUE)
ci(ct.bayes)
par(mfrow=c(2,2))
plot(ct)
plot(ct.bayes)
cc.bayes <- coloc.test(lm1,lm2, plots.extra=list(x=c("eta","theta"), y=c("lhood","lhood")),
                       bayes=TRUE, bayes.factor=list(c(-0.1,1), c(0.9,1.1)))
ci(cc.bayes)
## Bayesian approach, esp. when only p values are available
abf <- coloc.abf(list(beta=b1, varbeta=diag(v1), N=nrow(X1), sdY=sd(Y1), type="quant"),
                 list(beta=b2, varbeta=diag(v2), N=nrow(X2), sdY=sd(Y2), type="quant"))
abf
```

**garfield**

Again it can be installed with `biocLite("garfield")` and vignette be seen similarly to `coloc`.

> GWAS analysis of regulatory or functional information enrichment with LD correction. Briefly, it is a method that leverages GWAS findings with regulatory or 
> functional annotations (primarily from ENCODE and Roadmap epigenomics data) to find features relevant to a phenotype of interest. It performs greedy pruning of 
> GWAS SNPs (LD r2 > 0.1) and then annotates them based on functional information overlap. Next, it quantifies Fold Enrichment (FE) at various GWAS significance 
> cutoffs and assesses them by permutation testing, while matching for minor allele frequency, distance to nearest transcription start site and number of LD 
> proxies (r2 > 0.8).

**moloc**

moloc: multiple trait co-localization, available from https://github.com/clagiamba/moloc, can be installed with
```r
library(devtools)
install_github("clagiamba/moloc")
```

**rjags**

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

**sva**

The package contains function ```ComBat.R``` from https://www.bu.edu/jlab/wp-assets/ComBat/Download.html as described in the following paper.

> Johnson, WE, Rabinovic, A, and Li, C (2007). Adjusting batch effects in microarray expression data using Empirical Bayes methods. Biostatistics 8(1):118-127
```r
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
browseVignettes("sva")
```

**EBSeq**

The Bioconductor page is here, http://www.bioconductor.org/packages/devel/bioc/html/EBSeq.html.

## --- eQTL, epigenome-wide association study (EWAS). ---

**CpGassoc**

https://CRAN.R-project.org/package=CpGassoc

**missMethy**

http://bioconductor.org/packages/release/bioc/html/missMethyl.html

**QCGWAS**

It can be installed from CRAN. The sample is fairly easy to get going
```r
# 23-11-2018 JHZ

library(QCGWAS)
path <- "/home/jhz22/R/QCGWAS/data"
files <- file.path(path,dir(path))
load(files[1])
load(files[2])
head(gwa_sample,5)
head(header_translations,20)
write.table(gwa_sample,file="test",row.names=FALSE,quote=FALSE)
QCresults <- QC_GWAS("test",
		header_translations = header_translations,
		save_final_dataset = TRUE)
# allele-frequency threshold=0.05, HWEp=1e-4, call rate0.99, imputation quality=0.4
QCresults <- QC_GWAS("test",
		header_translations = header_translations,
		save_final_dataset = TRUE,
		HQfilter_FRQ = 0.05, HQfilter_HWE = 10^-4,
		HQfilter_cal = 0.99, HQfilter_imp = 0.4,
		NAfilter = TRUE)
# filters for the QQ-plot
QCresults <- QC_GWAS("test",
		header_translations = header_translations,
		save_final_dataset = TRUE,
		HQfilter_FRQ = 0.01, HQfilter_HWE = 10^-6,
		HQfilter_cal = 0.95, HQfilter_imp = 0.3,
		QQfilter_FRQ = c(NA, 0.01, 0.03, 0.05, 3),
		QQfilter_HWE = c(NA, 10^-6, 10^-4),
		QQfilter_cal = c(NA, 0.95, 0.98, 0.99),
		QQfilter_imp = c(NA, 0.3, 0.5, 0.7, 0.9),
		NAfilter = TRUE)
# HapMap allele reference -- but it does not work and should be
# https://ftp.hapmap.org/hapmap/frequencies/2010-08_phaseII+III/allele_freqs_chr2_CEU_r28_nr.b36_fwd.txt.gz
# based on https://bioinformatics.mdanderson.org/Software/VariantTools/mirror/annoDB/hapmap_CEU_freq.ann
# add options method="curl", extra="--insecure" to download.file
create_hapmap_reference(dir = ".",
		download_hapmap = TRUE,
		download_subset = "CEU",
		filename = "hapmap",
		save_txt = FALSE, save_rdata = TRUE)
# a new QC with HapMap
QCresults <- QC_GWAS("test",
		header_translations = header_translations,
		save_final_dataset = TRUE,
		HQfilter_FRQ = 0.01, HQfilter_HWE = 10^-6,
		HQfilter_cal = 0.95, HQfilter_imp = 0.3,
		QQfilter_FRQ = c(NA, 0.01, 0.03, 0.05, 3),
		QQfilter_HWE = c(NA, 10^-6, 10^-4),
		QQfilter_cal = c(NA, 0.95, 0.98, 0.99),
		QQfilter_imp = c(NA, 0.3, 0.5, 0.7, 0.9),
		NAfilter = TRUE,
		allele_ref_std = "hapmap.RData",
		allele_name_std = "HapMap",
		remove_mismatches = TRUE,
		check_ambiguous_alleles = FALSE)
# An alternative allele reference
QCresults <- QC_GWAS("test",
		header_translations = header_translations,
		save_final_dataset = TRUE,
		HQfilter_FRQ = 0.01, HQfilter_HWE = 10^-6,
		HQfilter_cal = 0.95, HQfilter_imp = 0.3,
		QQfilter_FRQ = c(NA, 0.01, 0.03, 0.05, 3),
		QQfilter_HWE = c(NA, 10^-6, 10^-4),
		QQfilter_cal = c(NA, 0.95, 0.98, 0.99),
		QQfilter_imp = c(NA, 0.3, 0.5, 0.7, 0.9),
		NAfilter = TRUE,
		allele_ref_std = "hapmap.RData",
		allele_name_std = "HapMap",
		remove_mismatches = TRUE,
		allele_ref_alt = NULL,
		allele_name_alt = "alternative",
		update_alt = TRUE,
		update_savename = "ref_alternative",
		update_as_rdata = TRUE)
# and QC with it
QCresults <- QC_GWAS("test",
		header_translations = header_translations,
		save_final_dataset = TRUE,
		HQfilter_FRQ = 0.01, HQfilter_HWE = 10^-6,
		HQfilter_cal = 0.95, HQfilter_imp = 0.3,
		QQfilter_FRQ = c(NA, 0.01, 0.03, 0.05, 3),
		QQfilter_HWE = c(NA, 10^-6, 10^-4),
		QQfilter_cal = c(NA, 0.95, 0.98, 0.99),
		QQfilter_imp = c(NA, 0.3, 0.5, 0.7, 0.9),
		NAfilter = TRUE,
		allele_ref_std = "hapmap.RData",
		allele_name_std = "HapMap",
		remove_mismatches = TRUE,
		allele_ref_alt = "ref_alternative.RData",
		allele_name_alt = "alternative",
		update_alt = TRUE,
		update_as_rdata = TRUE,
		backup_alt = TRUE)
# automatic loading
hapmap_ref <- read.table("hapmap_ref.txt",
		header = TRUE, as.is = TRUE)
alternative_ref <- read.table("alt_ref.txt",
		header = TRUE, as.is = TRUE)
QCresults <- QC_GWAS("test",
		header_translations = "headers.txt",
		out_header = "new_headers.txt",
		allele_ref_std = hapmap_ref,
		allele_ref_alt = alternative_ref,
		update_alt = TRUE,
		update_as_rdata = FALSE,
		update_savename = "alt_ref")
# automatic QC of multiple files
QC_series(
	data_files= c("data1.txt","data2.txt","data3.txt"),
	output_filenames = c("output1.txt","output2.txt","output3.txt"),
	dir_data = "preQC",
	dir_output = "postQC",
	dir_references = "QC_files",
	header_translations = header_translations,
	save_final_dataset = TRUE,
	HQfilter_FRQ = 0.01, HQfilter_HWE = 10^-6,
	HQfilter_cal = 0.95, HQfilter_imp = 0.3,
	QQfilter_FRQ = c(NA, 0.01, 0.03, 0.05, 3),
	QQfilter_HWE = c(NA, 10^-6, 10^-4),
	QQfilter_cal = c(NA, 0.95, 0.98, 0.99),
	QQfilter_imp = c(NA, 0.3, 0.5, 0.7, 0.9),
	NAfilter = TRUE,
	allele_ref_std = "ref_hapmap.RData",
	allele_name_std = "HapMap",
	remove_mismatches = TRUE,
	allele_ref_alt = "ref_alternative.RData",
	allele_name_alt = "alternative",
	update_alt = TRUE, update_as_rdata = TRUE, backup_alt = TRUE)
```
Note the changes required with the HapMap reference, i.e.,
```r
download.file(url = paste0("https://ftp.hapmap.org/hapmap/frequencies/2010-08_phaseII+III/", 
              "allele_freqs_chr", dn, "_", download_subset, 
              "_r28_nr.b36_fwd.txt.gz"), method = "curl", extra = "--insecure", 
              destfile = paste0(dir, "/allele_freqs_chr", dn, 
              "_", download_subset, "_r28_nr.b36_fwd.txt.gz"))
```
ftp://ftp.ncbi.nlm.nih.gov/hapmap/frequencies/2010-08_phaseII+III/ might also work.
