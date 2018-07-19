# Benchmarks and pipelines

## Agotron detection

The following is according to https://github.com/ncrnalab/agotron_detector as described in
> Hansen TB (2018). Detecting Agotrons in Ago CLIPseq Data. in Vang Ørom UA (ed) miRNA Biogenesis-Methods and Protocols, Chapter 17, 221-232. Springer.
```bash
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -zxvf chromFa.tar.gz
cat *.fa > hg19.fa
samtools faidx hg19.fa
bowtie2-build hg19.fa hg19
# GSE78059
for srr in 008/SRR3177718/SRR3177718 009/SRR3177719/SRR3177719 000/SRR3177720/SRR3177720 001/SRR3177721/SRR3177721 002/SRR3177722/SRR3177722 003/SRR3177723/SRR3177723
do
   wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR317/$srr.fastq.gz
done
trim_galore -A TCAGTCACTTCCAGC -length 18 *.fastq.gz
for i in *_trimmed.fq.gz
do
    echo $i
    bowtie2 -q --local -x hg19 -U $i | samtools sort - > $i.sort.bam    
    samtools index $i.sort.bam
    
done
python UCSC_intron_retriever.py | python analyzer.py -g hg19.fa | Rscript annotater.R
```
Note that it is easier to implement with ```prefetch``` as shown below.

## SNP discovery

The following reference discribes several pipelines for SNP discovery.

> Morin PA, Foote AD, Hill CM, Simon-Bouhet B, Lang AR, Louis M (2018). SNP Discovery from Single and Multiplex Genome Assemblies of Non-model Organisms, in Head SR, et al. (eds.), Next Generation Sequencing: Methods and Protocols, Chapter 9, 113-144, Springer.

whose scripts are available from https://github.com/PAMorin/SNPdiscovery/.

See also https://github.com/sanger-pathogens/snp-sites and the following references,

> Martin J, Schackwitz W, Lipzen A (2018). Genomic Sequence Variation Analysis by Resequencing, in de Vries RP, Tsang A, Grigoriev IV (ed) Fungal Genomics-Methods and Protocols, 2e, Chapter 18, 229-239, Springer.

> Raghavachari N, Garcia-Reyero N (eds.) (2018), Gene Expression Analysis-Methods and Protocols, Springer.

## GSNAP, MapSplice, RUM, STAR, RNA-seq pipeline

```bash
# gsnap
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-07-04.tar.gz
tar xfz gmap-gsnap-2018-07-04.tar.gz
cd gmap-2018-07-04
./configure
make
sudo make install
# mapsplice, the latest version from http://protocols.netlab.uky.edu/~zeng/MapSplice-v2.2.1.zip has compiling issue
sudo `which conda` install mapsplice
mapsplice.py
# rum
git clone https://github.com/itmat/rum
cd rum
perl Makefile.PL
make
sudo make install
# STAR
git clone https://github.com/alexdobin/STAR
cd STAR/source
make
```
See https://github.com/sanger-pathogens/Bio-RNASeq.

## sra-toolkit, tophat

These are very straightforward, e.g.,
```bash
prefetch -v SRR3534842
fastq-dump --split-files --gzip SRR3534842
```
the SRR3534842.sra from prefetch is actually at $HOME/ncbi/public/sra which is split into `SRR3534842_1.fastq.gz`, `SRR3534842_2.fastq.gz` at the current directory. See https://www.biostars.org/p/111040/. However, the location may not desirable since it may create a huge .vdi files with VirtualBox -- to get around we do this
```bash
cd $HOME
mkdir -p /home/jhz22/D/work/ncbi/public/sra
ln -sf /home/jhz22/D/work/ncbi
```
where D is actually a shared folder at Windows.

To run ```tophat```, see https://ccb.jhu.edu/software/tophat/tutorial.shtml
```bash
wget https://ccb.jhu.edu/software/tophat/downloads/test_data.tar.gz
tar xvfz test_data.tar.gz
cd test_data
tophat -r 20 test_ref reads_1.fq reads_2.fq
```

# Software

## ANGSD

ANGSD is a software for analyzing next generation sequencing data, http://www.popgen.dk/angsd/index.php/ANGSD. It is relatively straightforward with GitHub; after
```bash
git clone https://github.com/ANGSD/angsd
cd angsd
make
```
but the following change is needed on line 468 of `misc/msHOT2glf.c`: `tmppch` as in (tmppch=='\0') should be `*tmppch` as in `(*tmppch==''0')`, suggested by the compiler.

## bamtools, bcftools, bedops, bedtools, blast (ncbi-blast+), bowtie2, fastqc, fastx-toolkit, hmmer, hisat2, picard-tools, rsem, sambamba, samtools, seqtk, sra-toolkit, tophat, vcftools (Ubuntu archive)

Install with ```sudo apt install```.

See also https://github.com/lh3/seqtk.

# bowtie2

The project home is https://sourceforge.net/projects/bowtie-bio, whereby
```bash
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-source.zip
unzip bowtie2-2.3.4.1-linux-x86_64.zip
cd bowtie2-2.3.4.1-linux-x86_64/
```
The test is then self-contained,
```bash
export BT2_HOME=/home/jhz22/D/genetics/bowtie2-2.3.4.1-linux-x86_64

$BT2_HOME/bowtie2-build $BT2_HOME/example/reference/lambda_virus.fa lambda_virus
$BT2_HOME/bowtie2 -x lambda_virus -U $BT2_HOME/example/reads/reads_1.fq -S eg1.sam
$BT2_HOME/bowtie2 -x $BT2_HOME/example/index/lambda_virus -1 $BT2_HOME/example/reads/reads_1.fq -2 $BT2_HOME/example/reads/reads_2.fq -S eg2.sam

samtools view -bS eg2.sam > eg2.bam
samtools sort eg2.bam -o eg2.sorted.bam
samtools mpileup -uf $BT2_HOME/example/reference/lambda_virus.fa eg2.sorted.bam | bcftools view -Ov - > eg2.raw.bcf
bcftools view eg2.raw.bcf
```
Like samtools, etc. it is possible to involve `sudo apt install bowtie2`.

## cutadapt, TrimGalore

A prerequesite is to install cython.
```bash
git clone https://github.com/marcelm/cutadapt
cd cutadapt
sudo python setup.py install
git clone https://github.com/FelixKrueger/TrimGalore
```

## Exomiser

```bash
git clone https://github.com/exomiser/Exomiser
cd Exomiser
mvn package
```
See also https://github.com/exomiser/exomiser-demo.

## fastq-splitter

The scripts divides a large FASTQ file into a set of smaller equally sized files, http://kirill-kryukov.com/study/tools/fastq-splitter/.

## fastx_toolkit, RSEM

It is also available from https://github.com/agordon/fastx_toolkit along with https://github.com/agordon/libgtextutils, and do away with the notorious automake-1.14 problem associated with sources at http://hannonlab.cshl.edu/fastx_toolkit/download.html.

However, line 105 of ```src/fasta_formatter/fasta_formatter.cpp``` requires ```usage()``` followed by ```exit(0);``` as suggested in the `issue` section. More oever, usage() is a void function so its own `exit(0)` is unnecessary.

The GitHub pages for RSEM are https://github.com/deweylab/RSEM and https://deweylab.github.io/RSEM/. It is also recommended that the Bioconductor package EBSeq be installed.

## GATK

The source is available from https://github.com/broadinstitute/gatk/ but it is more convenient to use https://github.com/broadinstitute/gatk/releases/.
```bash
ln -s `pwd`/gatk $HOME/bin/gatk
gatk --help
gatk --list
```

## hisat2, sambamba, picard-tools, StringTie

Except StringTie, this is overlapped with ```apt install``` above,
```bash
brew tap brewsci/bio
brew tap brewsci/science
brew install hisat2
hisat2-build
brew install sambamba
brew install picard-tools
brew install stringtie
```
It could be useful with ``brew reinstall```. See
> Raghavachari N, Garcia-Reyero N (eds.) (2018), Gene Expression Analysis-Methods and Protocols, https://www.springer.com/us/book/9781493978335, Chapter 15, Springer.

Nevertheless it may be slower, e.g., tophat, compared to ```sudo apt install```.

## IGV

The download can be seeded from http://data.broadinstitute.org/igv, e.g., http://data.broadinstitute.org/igv/projects/downloads/2.4/IGV_2.4.10.zip.

Again the source code is from GitHub, https://github.com/igvteam/igv/. For developers, [igv.js](https://github.com/igvteam/igv.js) is very appealing.

## Jannovar

From the GitHub repository, it is seen to use `project object model` (POM), an XML representation of a Maven project held in a file named `pom.xml`. We therefore install `maven` first,
```bash
sudo apt install maven
```
The installation then proceeds as follows,
```bash
git clone https://github.com/charite/jannovar
cd jannovar
mvn package
```
Other tasks such as compile, test, etc. are also possible.

It is handy to use symbolic link, i.e.,
```bash
ln -s /home/jhz22/D/genetics/jannovar/jannovar-cli/target/jannovar-cli-0.24.jar $HOME/bin/Jannovar.jar
java -jar $HOME/bin/Jannovar.jar db-list
java -jar $HOME/bin/Jannovar.jar download -d hg19/refseq
```
We may need to set memory size, e.g., 
```bash
java -Xms2G -Xmx4G -jar $HOME/bin/Jannovar.jar
```

# pindel

The software can be obtained from https://github.com/genome/pindel.

After htslib is installed, the canonical instruction is to issue
```bash
git clone https://github.com/samtools/htslib
cd htslib
make
sudo make install
cd -
git clone https://github.com/genome/pindel
cd pindel
./INSTALL ../htslib
```
It is 'standard' to have complaints about pindel.cpp, bddate.cpp and genotyping.cpp,
for `abs()` rather than `fabs()` from the header file `cmath` have been used. The
issue goes away when `abs` is replaced with `fabs` and in the case of bddata.cpp,
it is also necessary to invoke the header, i.e.,
```cpp
#include <cmath>
```

## SnpEff, SnpSift, clinEff

It is straightforward with the compiled version from sourceforge, which also includes clinEff.
```bash
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core
cd snpEff
java -jar snpEff.jar databases
java -jar snpEff.jar download GRCh38.76
wget http://sourceforge.net/projects/snpeff/files/databases/test_cases.tgz
tar fvxz test_cases.tgz
```
lists all the databases and download a particular one. Later, the test files are also downloaded and extracted.

The following steps compile from source instead.
```bash
git clone https://github.com/pcingola/SnpEff.git
cd SnpEff
mvn package
mvn install
cd lib
# Antlr
mvn install:install-file \
	-Dfile=antlr-4.5.1-complete.jar \
	-DgroupId=org.antlr \
	-DartifactId=antlr \
	-Dversion=4.5.1 \
	-Dpackaging=jar

# BioJava core
mvn install:install-file \
	-Dfile=biojava3-core-3.0.7.jar \
	-DgroupId=org.biojava \
	-DartifactId=biojava3-core \
	-Dversion=3.0.7 \
	-Dpackaging=jar

# BioJava structure
mvn install:install-file \
	-Dfile=biojava3-structure-3.0.7.jar \
	-DgroupId=org.biojava \
	-DartifactId=biojava3-structure \
	-Dversion=3.0.7 \
	-Dpackaging=jar

cd -

# SnpSift
git clone https://github.com/pcingola/SnpSift.git
cd SnpSift
mvn package
mvn install
```
which gives `target/SnpEff-4.3.jar` and `target/SnpSift-4.3.jar`, respectively.

Note that `antlr4` is from GitHub, https://github.com/antlr/antlr4. See also https://github.com/sanger-pathogens/SnpEffWrapper.

# VarScan

Hosted at https://github.com/dkoboldt/varscan, the .jar files are ready to use with

```bash
git clone https://github.com/dkoboldt/varscan
```
or from the repository releases.

See http://varscan.sourceforge.net/ for further information.
