# Association analysis

## Data management

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
See [https://bitbucket.org/gavinband/bgen/overview](https://bitbucket.org/gavinband/bgen/overview).

## Single variant analysis

### eigensoft

The PCA software for genomewide data is available from [https://www.hsph.harvard.edu/alkes-price/software/](https://www.hsph.harvard.edu/alkes-price/software/) as well as Ubuntu.
```bash
sudo apt install eigensoft
```
The executables are eigenstrat, eigenstratQTL, smarteigenstrat, smartpca, pca, etc.

### GEMMA

To build from source, [https://github.com/genetics-statistics/GEMMA](https://github.com/genetics-statistics/GEMMA), the Makefile needs to change in places with OpenBLAS, /opt/OpenBLAS/.

### METAL

Note METAL aligns alleles according to the first file processed.

At least cmake 3.1 is required for the latest from GitHub, [https://github.com/statgen/METAL](https://github.com/statgen/METAL), 

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

As with distribution 2011-03-25, [http://csg.sph.umich.edu/abecasis/Metal/download/](http://csg.sph.umich.edu/abecasis/Metal/download/), options CUSTOMVARIABLE uses an output format of %g, leading to scientific notation of position, which
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
Another extension relates to heterogeneity analysis, e.g., I<sup>2</sup> > 30 we require at least three studies each attaining P <= 0.05. In this case,
we extend the direction field as in 
```c
        direction[marker] = z == 0.0 ? '0' : (z > 0.0 ? '+' : '-');
```
to
```c
        direction[marker] = z == 0.0 ? '0' : (z > 0.0 ? '+' : '-');
        direction[marker] = (fabs(z) * sqrt(w) < 1.959964) ? direction[marker] : (z > 0.0 ? 'p' : 'n');
```
for both ProcessFile() and ReProcessFile().

It is then relatively easy to filter on meta-analysis statistics, `awk -f metal.awk 4E.BP1-1.tbl`,
where `metal.awk` has the following lines,
```bash
{
   d3=$13; gsub(/?/,"",d3)
   if (length(d3) >= 3 && $18 >= 3500)
   if ($12 > -9.30103) print;
   else {
      if ($14 < 30) print;
      else {
        d3n=d3; d3p=d3;
        gsub(/+|-|p/,"",d3n); gsub(/+|-|n/,"",d3p);
        if (length(d3n) >= 3 || length(d3p) >= 3) print;
      }
   }
}
# R
# > log10(5e-10)
# [1] -9.30103
# head -1 METAL/4E.BP1-1.tbl | sed 's|\t|\n|g' | awk '{print "#" NR,$1}'
#1 Chromosome
#2 Position
#3 MarkerName
#4 Allele1
#5 Allele2
#6 Freq1
#7 FreqSE
#8 MinFreq
#9 MaxFreq
#10 Effect
#11 StdErr
#12 log(P)
#13 Direction
#14 HetISq
#15 HetChiSq
#16 HetDf
#17 logHetP
#18 N
```

### METASOFT and ForestPMPlot

Available from [http://genetics.cs.ucla.edu/meta/](http://genetics.cs.ucla.edu/meta/)
```bash
mkdir METASOFT
wget -qO- http://genetics.cs.ucla.edu/meta/repository/2.0.1/Metasoft.zip | \
unzip Metasoft.zip
# The results and log are in files `out` and `log`.
java -jar Metasoft.jar -input example.txt
java -jar Metasoft.jar -input example.txt -log example.log -output example.out
cd -

Somehow the website is down on 26/5/2022, and a copy is made available from [here](files/METASOFT).
```
One can also work with [ForestPMPLot](http://genetics.cs.ucla.edu/ForestPMPlot/repository/1.0.3/ForestPMPlot.zip) similarly.

---

### mtag

https://github.com/omeed-maghzian/mtag

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

### seqMeta

seqMeta: Meta-Analysis of Region-Based Tests of Rare DNA Variants, https://cran.r-project.org/web/packages/seqMeta/index.html.

### Quanto

[https://preventivemedicine.usc.edu/download-quanto/](https://preventivemedicine.usc.edu/download-quanto/)

### QUICKTEST

See [https://wp.unil.ch/sgg/quicktest/](https://wp.unil.ch/sgg/quicktest/) for the latest with support for bgen v1.2.


### SNPTEST

Available from [https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html#download](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html#download), e.g., 
```bash
wget http://www.well.ox.ac.uk/~gav/resources/snptest_v2.5.4-beta3_CentOS6.6_x86_64_static.tgz
tar xvfz nptest_v2.5.4-beta3_CentOS6.6_x86_64_static.tgz
```
It is possible that one would get error messages
```
!! Error in function: PerVariantComputationManager::get_phenotypes(), argument(s): phenotype_spec=IL.18R1___Q13478:P.
!! Quitting.
!! Error (HaltProgramWithReturnCode):
```
but they would go away with `-method em` for instance.

### SUGEN

Genetic Association Analysis Under Complex Survey Sampling, [https://github.com/dragontaoran/SUGEN](https://github.com/dragontaoran/SUGEN).

### swiss

Software to help identify overlap between association scan results and GWAS hit catalogs. 

```bash
sudo apt install libz-dev
pip install git+https://github.com/welchr/swiss.git@v1.0.0
```
One may try options such as --install_options="--prefix="". In case the $HOME directory does not have sufficient space, one can issue `swiss --download-data` on a system that does, then upload,
```bash
rsync -av --partial .local/share/swiss login.hpc.cam.ac.uk:$HOME/.local/share
```
or select particular files,
```bash
sync -av --partial .local/share/swiss/data/ld/1000g.phase3.hg38.EUR.shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz* \
     login.hpc.cam.ac.uk:$HOME/.local/share/swiss/data/ld
```
but this could be tricked with making $HOME/.local/share a symbolic pointing to a directory that can hold more than ~40GB data (reminscent fo VEP), then the `pip` command above can be called with the `--user` option.

To test, follow these,
```bash
git clone https://github.com/statgen/swiss
cd swiss/test
swiss --list-files
swiss --list-ld-sources
swiss --list-gwas-cats
swiss --list-gwas-traits --gwas-cat data/gwascat_ebi_GRCh37p13.tab

swiss --assoc data/top_hit_is_gwas.tab --variant-col EPACTS --pval-col PVAL \
      --dist-clump --clump-dist 250000 --clump-p 5e-08 --skip-gwas --out test

swiss --assoc data/test_hg19.gz --multi-assoc --trait SM --build hg19 \
      --ld-clump-source 1000G_2014-11_EUR --ld-gwas-source 1000G_2014-11_EUR --gwas-cat data/gwascat_ebi_GRCh37p13.tab \
      --ld-clump --clump-p 1e-10 --out test

swiss --assoc data/test_hg38.gz --gwas-cat data/gwascat_ebi_GRCh38p7.tab --variant-col VARIANT \
      --chrom-col CHROM --pos-col POS --trait BMI --build hg38 \
      --ld-clump-source 1000G_2014-11_EUR --ld-gwas-source 1000G_2014-11_EUR \
      --ld-clump --clump-p 5e-08 --out test
```
and consult the online documentation.

## Population struction

### fineSTRUCTURE

[https://people.maths.bris.ac.uk/~madjl/finestructure/index.html](https://people.maths.bris.ac.uk/~madjl/finestructure/index.html)

## HLA imputation

### HLA*IMP:02/03

Download source from [https://oxfordhla.well.ox.ac.uk/hla/tool/main](https://oxfordhla.well.ox.ac.uk/hla/tool/main)

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

There appears to have a web counterpart called [HLA*IMP:03](http://imp.science.unimelb.edu.au/hla/).

### Others

See [https://cambridge-ceu.github.io/csd3/applications/CookHLA.html](https://cambridge-ceu.github.io/csd3/applications/CookHLA.html) on SNP2HLA, HIBAG and CookHLA.

* [arcasHLA](https://github.com/RabadanLab/arcasHLA)
* [hla-genotyper](https://pypi.org/project/hla-genotyper/)
* [HLAgm](http://paed.hku.hk/genome/software/HLAgm_wes.zip)
* [HLAreporter](http://paed.hku.hk/genome/software/HLAreporter.zip)

Huang, Y. et al. HLAreporter: a tool for HLA typing from next generation sequencing data. Genome Medicine 7, 25 (2015).

---

## Finemapping

### CAVIAR/eCAVIAR/MsCAVIAR

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
cd caviar/CAVIAR-C++
make
cd -
CAVIAR -l CAVIAR-C++/sample_data/50_LD.txt -z CAVIAR-C++/sample_data/50_Z.txt -o 50
CAVIAR -l CAVIAR-C++/sample_data/DDB1.top100.sig.SNPs.ld -z CAVIAR-C++/sample_data/DDB1.top100.sig.SNPs.ZScores -o 100
eCAVIAR -l CAVIAR-C++/sample_data/GWAS.ADGC.MC.AD.IGAP.stage1.hg19.chr.11.121344805.121517613.CHRPOSREFALT.LD.ld \
        -l CAVIAR-C++/sample_data/eQTL.CARDIOGENICS.MC.AD.IGAP.stage1.hg19.chr.11.121344805.121517613.CHRPOSREFALT.LD.ld \
        -z CAVIAR-C++/sample_data/GWAS.MC.AD.IGAP.stage1.hg19.chr.11.121344805.121517613.CHRPOSREFALT.Z.txt \
        -z eQTL.CARDIOGENICS.MC.AD.IGAP.stage1.hg19.chr.11.121344805.121517613.CHRPOSREFALT.ILMN_1810712.NM_015313.1.ARHGEF12.Z.txt \
        -o 75
```
It may be necessary to alter Makefile to point to appropriate -I -L for lapack, for instance.

MsCAVIAR is available from [https://github.com/nlapier2/MsCAVIAR](https://github.com/nlapier2/MsCAVIAR).

The references are

**CAVIAR**

> Hormozdiari F, Kostem E, Kang EY, Pasaniuc B, Eskin E (2014). Identifying causal variants at loci with multiple signals of association. *Genetics* 198(2):497-508.

**eCAVIAR**

> Hormozdiari F, van de Bunt M, Segrè AV, Li X, Joo JWJ, Bilow M, Sul JH, Sankararaman S, Pasaniuc B, Eskin E (2016). Colocalization of GWAS and eQTL Signals Detects Target Genes. *Am J Hum Genet* 99(6):1245-1260.

Both are available from [https://github.com/fhormoz/caviar](https://github.com/fhormoz/caviar).

**MsCAVIAR**

LaPierre N, Taraszka K, Huang H, He R, Hormozdiari F,; Eskin E (2021): Additional simulations and a table of the real data results.. PLOS Genetics. Journal contribution. https://doi.org/10.1371/journal.pgen.1009733.s001.

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

### finemap

finemap 1.4 is available from [http://www.christianbenner.com/finemap_v1.4_x86_64.tgz](http://www.christianbenner.com/finemap_v1.4_x86_64.tgz).
```bash
finemap_v1.4_x86_64 --sss --in-files example/master
```
We could also keep the screen output with `finemap_v1.4_x86_64 --sss --in-files example/master 2>&1 | tee test.log`.

LDSTORE v2.0 is available from [http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz](http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz).

```bash
ldstore_v2.0_x86_64 --in-files example/master --write-bdose --bdose-version 1.1
ldstore_v2.0_x86_64 --in-files example/master --write-bcor --read-bdose
ldstore_v2.0_x86_64 --in-files example/master --bcor-to-text

```

### gchromVar

Cell type specific enrichments using finemapped variants and quantitative epigenetic data,
[https://github.com/caleblareau/gchromVAR](https://github.com/caleblareau/gchromVAR); see [https://github.com/caleblareau/singlecell_bloodtraits](https://github.com/caleblareau/singlecell_bloodtraits) for examples.

It requires chromVarmotifs from [https://github.com/GreenleafLab/chromVARmotifs](https://github.com/GreenleafLab/chromVARmotifs), which requires gcc 5.2.0.

### JAM

**Setup**

The package is available from [https://github.com/pjnewcombe/R2BGLiMS](https://github.com/pjnewcombe/R2BGLiMS).

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

## Functional annotation

### ANNOVAR

See [https://cambridge-ceu.github.io/csd3/applications/ANNOVAR.html](https://cambridge-ceu.github.io/csd3/applications/ANNOVAR.html).

### DAVID

[https://david.ncifcrf.gov/](https://david.ncifcrf.gov/)

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

### GARFIELD

Web: [https://www.ebi.ac.uk/birney-srv/GARFIELD/](https://www.ebi.ac.uk/birney-srv/GARFIELD/)
```bash
wget -qO- https://www.ebi.ac.uk/birney-srv/GARFIELD/package-v2/garfield-v2.tar.gz | \
tar xvfz
wget -qO- https://www.ebi.ac.uk/birney-srv/GARFIELD/package-v2/garfield-data.tar.gz | \
tar xvfz
cd garfield-v2
bash garfield
```
Additional details are described in the documentation, [https://www.ebi.ac.uk/birney-srv/GARFIELD/documentation-v2/GARFIELD-v2.pdf](https://www.ebi.ac.uk/birney-srv/GARFIELD/documentation-v2/GARFIELD-v2.pdf).

### gnomAD

Web site: [https://gnomad.broadinstitute.org/downloads](https://gnomad.broadinstitute.org/downloads).

To use [gsutil](https://cloud.google.com/storage/docs/gsutil), following these steps,
```bash
module load python/3.7
virtualenv py37
source py37/bin/activate
pip install gsutil
gsutil ls gs://gnomad-public/release
gsutil cp gs://gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz .
gsutil cp gs://gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz .
```
UCSC has a description [here](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&c=chrX&g=gnomadPLI).

### HaploReg

Web: [https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php](https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php) (data files, [https://pubs.broadinstitute.org/mammals/haploreg/data/](https://pubs.broadinstitute.org/mammals/haploreg/data/))

> HaploReg is a tool for exploring annotations of the noncoding genome at variants on haplotype blocks, such as candidate regulatory SNPs at disease-associated loci. Using LD information from the 1000 Genomes Project, linked SNPs and small indels can be visualized along with chromatin state and protein binding annotation from the Roadmap Epigenomics and ENCODE projects, sequence conservation across mammals, the effect of SNPs on regulatory motifs, and the effect of SNPs on expression from eQTL studies. HaploReg is designed for researchers developing mechanistic hypotheses of the impact of non-coding variants on clinical phenotypes and normal variation.

### morpheus

See [https://software.broadinstitute.org/morpheus/](https://software.broadinstitute.org/morpheus/).

### PolyPhen and PolyPhen-2

See [http://genetics.bwh.harvard.edu/pph/](http://genetics.bwh.harvard.edu/pph/) and [http://genetics.bwh.harvard.edu/pph2/](http://genetics.bwh.harvard.edu/pph2/).

## UNPHASED

See [https://sites.google.com/site/fdudbridge/software/unphased-3-1](https://sites.google.com/site/fdudbridge/software/unphased-3-1).

### VEP

The description is available from [http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html](http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html).
```bash
git clone https://github.com/Ensembl/ensembl-vep.git

cd ensembl-vep
git pull
git checkout release/92
perl INSTALL.pl
```
The last line requires modules DBI, Build as described in the [LANGUAGES](https://jinghuazhao.github.io/Computational-Statistics/LANGUAGES/) section of [Computational-Statistics](https://jinghuazhao.github.io/Computational-Statistics/).

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

Recent notes on ANNOVAR and VEP are available from here, [https://cambridge-ceu.github.io/csd3/applications/VEP.html](https://cambridge-ceu.github.io/csd3/applications/VEP.html).

### R-packages

See R-packages section.

---

## Pathway analysis

### DEPICT

See the [GIANT+Biobank BMI analysis](https://jinghuazhao.github.io/Omics-analysis/BMI/).

**Installation and documentation example**

* The official site, [https://data.broadinstitute.org/mpg/depict/documentation.html](https://data.broadinstitute.org/mpg/depict/documentation.html) has links on [DEPICT_v1_rel194.tar.gz](https://data.broadinstitute.org/mpg/depict/depict_download/bundles/DEPICT_v1_rel194.tar.gz),
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

### ToppGene

A portal for gene list enrichment analysis and candidate gene prioritization based on functional annotations and protein interactions network

See [https://toppgene.cchmc.org/](https://toppgene.cchmc.org/).

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

## Mendelian randomisation

Bias and Type 1 error rate for Mendelian randomization with sample overlap, https://sb452.shinyapps.io/overlap/

Power calculation. [https://shiny.cnsgenomics.com/mRnd/](https://shiny.cnsgenomics.com/mRnd/).

### LCV

Software in Matlab and R for Latent Causal Variable model inferring genetically causal relationships using GWAS data.

[https://github.com/lukejoconnor/LCV](https://github.com/lukejoconnor/LCV)

See R-packages section.

## Polygenic modeling

### GCTA

It is possible to use dosage, a Bash function is as follows,
```bash
function INTERVAL_dosage()
{
  if [ ! -f nodup/${pr}.gen.gz ]; then qctool -g nodup/${pr}.bgen -og nodup/${pr}.gen.gz; fi

  gunzip -c nodup/${pr}.gen.gz | \
  awk -v sample=${sample} '
  {
     N=(NF-6)/3
     for(i=1;i<=N;i++) dosage[NR,i]=$((i-1)*3+8)+2*$((i-1)*3+9)
  } END {
     i=0;
     while (getline gf < sample) {
       split(gf,a);
       i++
       id[i]=a[1]
     }
     close(sample)
     for(i=1;i<=N;i++)
     {
       printf id[i+2] " ML_DOSE"; for(j=1;j<=NR;j++) printf " " dosage[j,i]; printf "\n"
     }
  }' | \
  gzip -f > nodup/${pr}.dosage.gz
  (
    echo SNP Al1 Al2 Freq1 MAF Quality Rsq
    qctool -g ${pr}.bgen -snp-stats -osnp - | \
    sed '1,9d' | \
    cut -f2,5,6,12,14,17,18 | \
    sed 's/\t/ /g;s/NA/0/g'
  ) | \
  grep -v Completed | \
  gzip -f > nodup/${pr}.info.gz
}

gcta-1.9 --dosage-mach-gz nodup/$pr.dosage.gz nodup/$pr.info.gz --make-grm-bin --out nodup/${pr}
```
where `pr' is input file root and `sample` is the associate sample file from which sample IDs are extracted.

When the imputed genotyeps are MaCH-based, it is possible to use `DosageConverter`.
```bash
## assuming you use hpc-work/ with a subdirectory called bin/
cd /rds/user/$USER/hpc-work/
git clone https://github.com/Santy-8128/DosageConvertor
cd DosageConverter
pip install cget --user
module load cmake-3.8.1-gcc-4.8.5-zz55m7x
./install.sh
cd /rds/user/$USER/hpc-work/bin/
ln -s /rds/user/$USER/hpc-work/DosageConvertor/release-build/DosageConvertor
## testing
DosageConvertor  --vcfDose  test/TestDataImputedVCF.dose.vcf.gz \
                 --info     test/TestDataImputedVCF.info \
                 --prefix   test \
                 --type     mach

gunzip -c test.mach.dose.gz | wc -l

DosageConvertor  --vcfDose  test/TestDataImputedVCF.dose.vcf.gz \
                 --info     test/TestDataImputedVCF.info \
                 --prefix   test \
                 --type     plink

gunzip -c test.plink.dosage.gz | wc -l
```
so the MaCH dosage file is individual x genotype whereas PLINK dosage file is genotype x individual.

### HESS

HESS (Heritability Estimation from Summary Statistics) is now available from [https://github.com/huwenboshi/hess](https://github.com/huwenboshi/hess) and has a web page at

[https://huwenboshi.github.io/hess-0.5/#hess](https://huwenboshi.github.io/hess-0.5/#hess)

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
python hess.py --prefix step1 --reinflate-lambda-gc 1 --tot-hsqg 0.8 0.2 --out step2
```
where snp150.txt from UCSC is described at the SUMSTATS repository, [https://github.com/jinghuazhao/SUMSTATS](https://github.com/jinghuazhao/SUMSTATS).

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

### MiXeR

[https://github.com/precimed/mixer](https://github.com/precimed/mixer)

### PLINK2

Both [PLINK 1.90 beta](https://www.cog-genomics.org/plink2/) and [PLINK 2.00 alpha](https://www.cog-genomics.org/plink/2.0/) have issue with .grm.bin.N which is shorter than expected for GCTA. The problem is insidious but would prevent chromosome-specific GRMs to be combined.

Nevertheless there is no such problem with its --make-grm-list which allows for the possibility to use --mgrm-list option to combine chromosome-specific GRMs.

Note also the way to use individual's IDs in PLINK2.

## PRrice-2

Source [https://github.com/choishingwan/PRSice](https://github.com/choishingwan/PRSice)

and documentation, [https://choishingwan.github.io/PRSice/](https://choishingwan.github.io/PRSice/) (scripts [pgs.sh](files/pgs.sh)).

see also [https://github.com/pgormley/polygenic-risk-scores](https://github.com/pgormley/polygenic-risk-scores).

### R-packages

See R-packages section.

---

## PheWAS

### CLARITE

See [https://github.com/HallLab](https://github.com/HallLab) for the R and Python packages, both linking RNHANES.

Lucas AM, et al. CLARITE facilitates the quality control and analysis process for EWAS of metabolic-related traits. *Fron Genet*.

See also wiki resources section of Omics-analysis as well as implementations in R-packages section.

---

## Transcriptome-wide association analysis (TWAS)

### IronThrone-GoT

[https://github.com/landau-lab/IronThrone-GoT](https://github.com/landau-lab/IronThrone-GoT)

Nam AS, et al. (2019). Somatic mutations and cell identity linked by Genotyping of Transcriptomes. *Nature" 571: 355–360.

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

This section follows [http://gusevlab.org/projects/fusion/](http://gusevlab.org/projects/fusion/) and is more compact. To install we do,
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

We could simplify the recent GEUV example script as follows,

```bash
#!/usr/bin/bash

export FUSION=${HPC_WORK}/fusion_twas
for d in work GEUV; do if [ ! -d ${FUSION}/${d} ]; then mkdir ${FUSION}/${d}; fi; done
cd ${FUSION}/work
rm *
ln -sf . output

export LDREF=/rds/user/jhz22/hpc-work/fusion_twas/LDREF
export PRE_GEXP=${HPC_WORK}/fusion_twas/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt
gunzip -c $PRE_GEXP.gz |  awk '! ($3 ~ "X") && NR >1' | while read PARAM;
do
  export CHR=$(echo $PARAM | awk '{ print $3 }')
  export P0=$(echo $PARAM | awk '{ print $4 - 0.5e6 }')
  export P1=$(echo $PARAM | awk '{ print $4 + 0.5e6 }')
  export OUT=$(echo $PARAM | awk '{ print $1 }')
  echo $CHR $P0 $P1 $OUT
  echo $PARAM | tr ' ' '\n' | tail -n+5 | paste <(gunzip -c $PRE_GEXP.gz | head -n1 | tr '\t' '\n' | tail -n+5 | awk '{ print $1,$1 }') - > $OUT.pheno
  plink2 --bfile $LDREF/1000G.EUR.$CHR --pheno $OUT.pheno --make-bed --out $OUT --chr $CHR --from-bp $P0 --to-bp $P1 > /dev/null
  Rscript ${FUSION}/FUSION.compute_weights.R --bfile $OUT --tmp $OUT.tmp --out $FUSION/GEUV/$OUT \
          --save_hsq --hsq_p 0.1 --models blup,lasso,top1,enet --verbose 2
done

# https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/analysis_results/
```
Note that `gcta_nr_robust`, `plink2` and `gemma` are already in the searching path and we tricked GEMMA with `ln -sf . output` with the current working directory.

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

### enloc

Available from [https://github.com/xqwen/integrative](https://github.com/xqwen/integrative). More recent version is [fastenloc](https://github.com/xqwen/fastenloc); also related are [dap](https://github.com/xqwen/dap/) and [torus](https://github.com/xqwen/torus/).

For instance torus can be installed as follows,
```bash 
git clone https://github.com/xqwen/torus/
cd torus
module load gsl
module load boost/1.49.0-gcc4.9.1
module load zlib
cd src
make
make static
mv torus torus.static ${HPC_work}/bin
cd -
```
and for the documentaion example on height, we have
```bash
torus -d Height.torus.zval.gz --load_zval -dump_pip Height.gwas.pip
gzip Height.gwas.pip

fastenloc -eqtl gtex_v8.eqtl_annot.vcf.gz -gwas Height.gwas.pip.gz -prefix Height
sort -grk6 Height.enloc.sig.out
```
Note that these use hg38 references provided, and it is possible to generate the hg19 counterparts via script [gtex_v8_hg19.sh](files/gtex_v8_hg19.sh).

### FINEMAP colocalization pipeline

This is available from [https://bitbucket.org/mgloud/production_coloc_pipeline](https://bitbucket.org/mgloud/production_coloc_pipeline), note also https://github.com/boxiangliu/locuscomparer.

### GWAS-PW

Available from [https://github.com/joepickrell/gwas-pw](https://github.com/joepickrell/gwas-pw). The installation is straightforward after boost library is available.

We can use the following code for the documentation example,
```bash
gwas-pw -i example_data/aam_height_example.gz -bed ${HPC_WORK}/ldetect/ldetect-data/EUR.bed -phenos AAM HEIGHT -o example_data/aam_height
```
where EUR.bed contains the information for approximately independent LD blocks.

### INFERNO / SparkINFERNO

Short for (INFERring the molecular mechanisms of NOncoding genetic variants, it is available from [https://bitbucket.org/wanglab-upenn/INFERNO](https://bitbucket.org/wanglab-upenn/INFERNO)
and it also has a web interface, [http://inferno.lisanwanglab.org/index.php](http://inferno.lisanwanglab.org/index.php).

SparkINFERNO is described here, [https://bitbucket.org/wanglab-upenn/sparkinferno/](https://bitbucket.org/wanglab-upenn/sparkinferno/).

Web [http://inferno.lisanwanglab.org/index.php](http://inferno.lisanwanglab.org/index.php).

### R packages

## biMM

It is a software for bivariate lineax mixed model (LMM),

[https://www.mv.helsinki.fi/home/mjxpirin/download.html](https://www.mv.helsinki.fi/home/mjxpirin/download.html)

pSI, available from CRAN and [http://genetics.wustl.edu/jdlab/psi_package/](http://genetics.wustl.edu/jdlab/psi_package/) with supplementary data [http://genetics.wustl.edu/jdlab/files/2014/01/pSI.data_1.0.tar_.gz](http://genetics.wustl.edu/jdlab/files/2014/01/pSI.data_1.0.tar_.gz).

## Epigenomics

**Avocado**

Project page, [https://noble.gs.washington.edu/proj/avocado/](https://noble.gs.washington.edu/proj/avocado/) and Software, [https://bitbucket.org/noblelab/avocado/src/master/](https://bitbucket.org/noblelab/avocado/src/master/)

To install
```bash
python setup.py install --prefix=/rds/user/jhz22/hpc-work
# insert this line into .bashrc
export PYTHONPATH=$PYTHONPATH:/rds/user/jhz22/hpc-work/lib/python2.7/site-packages/
```

**combined-pvalues**

A library to combine, analyze, group and correct p-values in BED files. Unique tools involve correction for spatial autocorrelation. This is useful for ChIP-Seq probes and Tiling arrays, or any data with spatial correlation.

Software, [https://github.com/brentp/combined-pvalues](https://github.com/brentp/combined-pvalues)

Pedersen BS, Schwartz DA, Yang IV, Kechris KJ. Comb-p: software for combining, analyzing, grouping and correcting spatially correlated P-values
*Bioinformatics* 28(22):2986–2988, https://doi.org/10.1093/bioinformatics/bts545

**ChromImpute**

Website, [http://www.biolchem.ucla.edu/labs/ernst/ChromImpute/](http://www.biolchem.ucla.edu/labs/ernst/ChromImpute/) and Source, [https://github.com/jernst98/ChromImpute](https://github.com/jernst98/ChromImpute)

Ernst J, Kellis M. Large-scale imputation of epigenomic datasets for systematic annotation of diverse human tissues. Nature Biotechnology, 33:364-376, 2015

**EpiAlign**

Software, [https://github.com/zzz3639/EpiAlign](https://github.com/zzz3639/EpiAlign)

Web, [http://shiny.stat.ucla.edu:3838/EpiAlign/](http://shiny.stat.ucla.edu:3838/EpiAlign/)

Ge X, Zhang H, Xie L, Li WV, Kwon SB, Li JJ
EpiAlign: an alignment-based bioinformatic tool for comparing chromatin state sequences.
Nucleic Acids Res. 2019 Apr 24. doi: 10.1093/nar/gkz287.

---

## R-packages

See [https://github.com/jinghuazhao/Computational-Statistics](https://github.com/jinghuazhao/Computational-Statistics) for general information.

**Bioconductor**

This includes Biobase, BSGenome, edgeR, limma, Rsubread, STRINGdb.

**CRAN**

This includes DCGL.

**GSMR**

It refers to Generalised Summary-data-based Mendelian Randomisation, http://cnsgenomics.com/software/gsmr/, available both as part of GCTA and R package.
```r
install.packages("http://cnsgenomics.com/software/gsmr/static/gsmr_1.0.6.tar.gz",repos=NULL,type="source")
```
with test data, [http://cnsgenomics.com/software/gsmr/static/test_data.zip](http://cnsgenomics.com/software/gsmr/static/test_data.zip).

Script for the documentation example is tallied here as [gsmr_example.R](files/gsmr_example.R).

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

### BMI and T2D
```r
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

See [https://github.com/PheWAS/](https://github.com/PheWAS/)

**coloc**

It requires `snpStats` that can be installed with biocLite().

Now it has a comprehensive site, [https://chr1swallace.github.io/coloc/](https://chr1swallace.github.io/coloc/)

This is the original example from the documentation,
```r
set.seed(1)
X1 <- matrix(rbinom(1200,1,0.4),ncol=2)
X2 <- matrix(rbinom(1000,1,0.6),ncol=2)
colnames(X1) <- colnames(X2) <- c("f1","f2")
library(dplyr)
X1 <- mutate(as.data.frame(X1),x0=1) %>% select(x0,f1,f2)
X2 <- mutate(as.data.frame(X2),x0=1) %>% select(x0,f1,f2)
Y1 <- rnorm(600,apply(X1,1,sum),2)
Y2 <- rnorm(500,2*apply(X2,1,sum),5)
summary(lm1 <- lm(Y1~f1+f2,data=as.data.frame(X1)))
summary(lm2 <- lm(Y2~f1+f2,data=as.data.frame(X2)))
b1 <- coef(lm1)
b2 <- coef(lm2)
v1 <- vcov(lm1)
v2 <- vcov(lm2)
require(coloc)
## Bayesian approach, esp. when only p values are available
abf <- coloc.abf(list(beta=b1, varbeta=diag(v1), N=nrow(X1), sdY=sd(Y1), type="quant"),
                 list(beta=b2, varbeta=diag(v2), N=nrow(X2), sdY=sd(Y2), type="quant"))
abf
# sdY
cat("sd(Y)=",sd(Y1),"==> Estimates:",sqrt(diag(var(X1)*b1^2+var(X1)*v1*nrow(X1))),"\n")
for(k in 1:3) cat("Based on b",k," sd(Y1) = ",sqrt(var(X1[,k])*(b1[k]^2+nrow(X1)*v1[k,k])),"\n",sep="")
cat("sd(Y)=",sd(Y2),"==> Estimates:",sqrt(diag(var(X2)*b2^2+var(X2)*v2*nrow(X2))),"\n")
for(k in 1:3) cat("Based on b",k," sd(Y2) = ",sqrt(var(X2[,k])*(b2[k]^2+nrow(X2)*v2[k,k])),"\n",sep="")
legacy <- function()
## intuitive test for proportionality from https://cran.r-project.org/web/packages/coloc/vignettes/vignette.html
{
# coloc, large (>0.05) p.value.chisquare indicates traits are compatible with colocalisation
  ct <- coloc.test(lm1,lm2, plots.extra=list(x=c("eta","theta"), y=c("lhood","lhood")))
  summary(ct)
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
}
```
where we have illustrated how to obtain sd(Y) whose outputs are as follows,
```
> # sdY
> cat("sd(Y)=",sd(Y1),"==> Estimates:",sqrt(diag(var(X1)*b1^2+var(X1)*v1*nrow(X1))),"\n")
sd(Y)= 2.087048 ==> Estimates: 0 2.070751 2.041013
> for(k in 1:3) cat("Based on b",k," sd(Y1) = ",sqrt(var(X1[,k])*(b1[k]^2+nrow(X1)*v1[k,k])),"\n",sep="")
Based on b1 sd(Y1) = 0
Based on b2 sd(Y1) = 2.070751
Based on b3 sd(Y1) = 2.041013
> cat("sd(Y)=",sd(Y2),"==> Estimates:",sqrt(diag(var(X2)*b2^2+var(X2)*v2*nrow(X2))),"\n")
sd(Y)= 5.694518 ==> Estimates: 0 5.592806 5.59157
> for(k in 1:3) cat("Based on b",k," sd(Y2) = ",sqrt(var(X2[,k])*(b2[k]^2+nrow(X2)*v2[k,k])),"\n",sep="")
Based on b1 sd(Y2) = 0
Based on b2 sd(Y2) = 5.592806
Based on b3 sd(Y2) = 5.59157
```
In fact, we usually only works with single regression coefficients.

Developmental version of the package is available as follows,
```r
if(!require("remotes"))
   install.packages("remotes")
remotes::install_github("chr1swallace/coloc")
```

**garfield**

Web site: [https://www.ebi.ac.uk/birney-srv/GARFIELD/](https://www.ebi.ac.uk/birney-srv/GARFIELD/)

Again it can be installed with `biocLite("garfield")` and vignette be seen similarly to `coloc`.

> GWAS analysis of regulatory or functional information enrichment with LD correction. Briefly, it is a method that leverages GWAS findings with regulatory or 
> functional annotations (primarily from ENCODE and Roadmap epigenomics data) to find features relevant to a phenotype of interest. It performs greedy pruning of 
> GWAS SNPs (LD r2 > 0.1) and then annotates them based on functional information overlap. Next, it quantifies Fold Enrichment (FE) at various GWAS significance 
> cutoffs and assesses them by permutation testing, while matching for minor allele frequency, distance to nearest transcription start site and number of LD 
> proxies (r2 > 0.8).

The documentation example is run as follows,
```r
     garfield.run("tmp", data.dir=system.file("extdata",package = "garfield"),
         trait="trait",run.option = "prep", chrs = c(22),
         exclude = c(895, 975, 976, 977, 978, 979, 98))

     garfield.run("tmp", data.dir=system.file("extdata",package = "garfield"),
         run.option = "perm", nperm = 1000, thresh = c(0.001, 1e-04, 1e-05),
         pt_thresh = c(1e-04, 1e-05), maf.bins = 2, tags.bins = 3, tss.bins = 3,
         prep.file = "tmp.prep", optim_mode = TRUE, minit = 100, thresh_perm = 0.05)

     if (file.exists("tmp.perm")){
         perm = read.table("tmp.perm", header=TRUE)
         head(perm)
     } else { print("Error: tmp.perm does not exist!") }
```

We have the Crohn's disease example,
```r
# download data and decompress
system("wget https://www.ebi.ac.uk/birney-srv/GARFIELD/package/garfield-data.tar.gz")
system("tar -zxvf garfield-data.tar.gz")

# if downloaded in current working directory use the following to execute
# garfield, otherwise please change data.dir location
garfield.run("cd-meta.output", data.dir="garfield-data", trait="cd-meta",
             run.option = "prep", chrs = c(1:22), exclude = c(895, 975, 976, 977, 978,
             979, 980))
#
garfield.run("cd-meta.output", data.dir="garfield-data", run.option = "perm",
             nperm = 100000, thresh = c(0.1,0.01,0.001, 1e-04, 1e-05, 1e-06, 1e-07, 1e-08),
             pt_thresh = c(1e-05, 1e-06, 1e-07, 1e-08), maf.bins = 5, tags.bins = 5,
             tss.bins = 5, prep.file = "cd-meta.output.prep", optim_mode = TRUE,
             minit = 100, thresh_perm = 0.0001)
#
garfield.plot("cd-meta.output.perm", num_perm = 100000,
              output_prefix = "cd-meta.output", plot_title = "Crohn's Disease",
              filter = 10, tr = -log10(0.05/498))
```

**GenomicSEM**

[GenomicSEM](https://github.com/MichelNivard/GenomicSEM) fits structural equation models based on the summary statistics obtained from genome wide association studies (GWAS).

**gtx**

The version at [https://github.com/tobyjohnson/gtx](https://github.com/tobyjohnson/gtx) has more updates to its counterpart at CRAN.

**HIBAG**

Currently it is archived at CRAN but can be downloaded from GitHub, [https://github.com/cran/HIBAG](https://github.com/cran/HIBAG)
```R
devtools::install_github("cran/HIBAG")
```
However it is now available from Bioconductor.

**hyprcoloc**

It is a package for hypothesis prioritisation multi-trait colocalization, available from [https://github.com/jrs95/hyprcoloc](https://github.com/jrs95/hyprcoloc).
```bash
R --no-save -q <<END
  library(hyprcoloc)
  hyprcoloc_test <- function()
  {
  # Regression coefficients and standard errors from ten GWAS studies (Traits 1-5, 6-8 & 9-10 colocalize)
    betas <- hyprcoloc::test.betas
    head(betas)
    ses <- hyprcoloc::test.ses
    head(ses)
  # Trait names and SNP IDs
    traits <- paste0("T", 1:10)
    rsid <- rownames(betas)
  # Colocalisation analyses
    hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
  }
  hyprcoloc_test()
END
```
**meta**

The following code, courtesy of the package developer, generates three forest plots,
```R
library(meta)

## Generic inverse-variance meta-analysis
## (first two arguments: treatment estimate and its standard error)
##
m1 <- metagen(1:10, rep(0.1, 10), sm = "MD", studlab = LETTERS[1:10])

## Use update.meta() to re-run meta-analysis with additional argument
##
m1.subset <- update(m1, subset = 1:5)
##
m1.exclude <- update(m1, exclude = 6:10)

pdf("forest1-all.pdf", width = 8.75, height = 4)
forest(m1, colgap.forest.left = "1cm")
grid::grid.text("All studies", 0.5, 0.94, gp = grid::gpar(cex = 1.5))

forest(m1.subset, colgap.forest.left = "1cm")
grid::grid.text("Subset of studies", 0.5, 0.9, gp = grid::gpar(cex = 1.5))

forest(m1.exclude, colgap.forest.left = "1cm")
grid::grid.text("Exclude studies", 0.5, 0.94, gp = grid::gpar(cex = 1.5))
dev.off()
```

**moloc**

moloc: multiple trait co-localization, available from [https://github.com/clagiamba/moloc](https://github.com/clagiamba/moloc), can be installed with
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

The package contains function ```ComBat.R``` from [https://www.bu.edu/jlab/wp-assets/ComBat/Download.html](https://www.bu.edu/jlab/wp-assets/ComBat/Download.html) as described in the following paper.

> Johnson, WE, Rabinovic, A, and Li, C (2007). Adjusting batch effects in microarray expression data using Empirical Bayes methods. Biostatistics 8(1):118-127
```r
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
browseVignettes("sva")
```

**EBSeq**

The Bioconductor page is here, [http://www.bioconductor.org/packages/devel/bioc/html/EBSeq.html](http://www.bioconductor.org/packages/devel/bioc/html/EBSeq.html).

## --- eQTL, epigenome-wide association study (EWAS) ---

**CpGassoc**

[https://CRAN.R-project.org/package=CpGassoc](https://CRAN.R-project.org/package=CpGassoc)

**fastQTL**

[http://fastqtl.sourceforge.net/files/FastQTL-2.165.linux.tgz](http://fastqtl.sourceforge.net/files/FastQTL-2.165.linux.tgz)

**missMethy**

[http://bioconductor.org/packages/release/bioc/html/missMethyl.html](http://bioconductor.org/packages/release/bioc/html/missMethyl.html)

**MultiABEL**

Multi-Trait Genome-Wide Association Analysis

[https://github.com/xiashen/MultiABEL](https://github.com/xiashen/MultiABEL)

**OmnibusFisher**

The p-values of SNPs, RNA expressions and DNA methylations are calculated by kernel machine (KM) regression. The correlation between different omics data are taken into account. This method can be applied to either samples with all three types of omics data or samples with two types.

[https://CRAN.R-project.org/package=OmnibusFisher](https://CRAN.R-project.org/package=OmnibusFisher)


**PredictABEL**

See [https://CRAN.R-project.org/package=PredictABEL](https://CRAN.R-project.org/package=PredictABEL).

Kundu S, et al. (2014). Estimating the predictive ability of genetic risk models in simulated data based on published results from genome-wide association studies. *Front Genet* 2014, 5: 179, [https://doi.org/10.3389/fgene.2014.00179](https://doi.org/10.3389/fgene.2014.00179).

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

**SAIGE**

The address of GitHub repository is here, [https://github.com/weizhouUMICH](https://github.com/weizhouUMICH)

Information including installation is described here, [https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE).

Locally, we therefore check for the desired gcc, cmake and boost and proceed as follows,
```bash
module avail gcc
module avail cmake
module avail boost
# gcc > 5.5 and cmake > 3.8.1
module load gcc-6.1.0-gcc-4.8.5-jusvegv cmake-3.8.1-gcc-4.8.5-zz55m7
# boost 1.58.0 and R 3.6.0 are described in Computationl_Statistics repository.
export LD_LIBRARY_PATH=/rds-d4/user/jhz22/hpc-work/boost_1_58_0/stage/lib:$LD_LIBRARY_PATH
# we actually use the binary disrtibution directly
wget https://github.com/weizhouUMICH/SAIGE/archive/v0.35.8.2.tar.gz
tar tvfz v0.35.8.2.tar.gz
cd SAIGE-0.35.8.2
tar xvfz SAIGE_0.35.8.2_R_x86_64-pc-linux-gnu.tar.gz
mv SAIGE /rds-d4/user/jhz22/hpc-work/R
R --no-save <<END
  library(SAIGE)
END
```
Since the required packages Rcpp and RcppParallel are relatively easy to deal with, with which we then simply load the packague as usual.

A recent description is given here, [https://cambridge-ceu.github.io/csd3/applications/SAIGE.html](https://cambridge-ceu.github.io/csd3/applications/SAIGE.html).

**ensemblVEP**

There is no particular difficulty, simply use `BiocManager::install("ensemblVEP")`.
