# Linuxbrew

Follow http://linuxbrew.sh/ and possibly https://docs.brew.sh
```bash
sudo apt-get install build-essential
sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
echo 'export PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"' >>~/.profile
echo 'export MANPATH="/home/linuxbrew/.linuxbrew/share/man:$MANPATH"' >>~/.profile
echo 'export INFOPATH="/home/linuxbrew/.linuxbrew/share/info:$INFOPATH"' >>~/.profile
PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"
```

# hisat2, sambamba, picard-tools, StringTie

```bash
brew tap brewsci/bio
brew tap brewsci/science
brew install hisat2
hisat2-build
brew install sambamba
brew install picard-tools
brew install stringtie
```
See https://www.springer.com/us/book/9781493978335, Chapter 15.

# tensorflow

The tensorflow repository is here, https://github.com/tensorflow/tensorflow, and it is relatively easy to install via pip,
```bash
pip install tensorflow
python <<END
import tensorflow as tf
hello = tf.constant('Hello, TensorFlow!')
sess = tf.Session()
print(sess.run(hello))
END
```
Follow https://github.com/aymericdamien/TensorFlow-Examples for readily adaptible examples.

Also

https://github.com/apress/pro-deep-learning-w-tensorflow

# sra-toolkit, samtools, bamtools, bcftools vcftools

These are very straightforward under Ubuntu, e.g.,
```bash
sudo apt install sra-toolkit
prefetch -v SRR3534842
fastq-dump --split-files --gzip SRR3534842
```
the SRR3534842.sra from prefetch is actually at $HOME/ncbi/public/sra which is split into
`SRR3534842_1.fastq.gz`, `SRR3534842_2.fastq.gz` at the current directory.

See https://www.biostars.org/p/111040/.

# cuadapt, TrimGalore

A prerequesite is to install cython.
```bash
git clone https://github.com/marcelm/cutadapt
cd cutadapt
sudo python setup.py install
git clone https://github.com/FelixKrueger/TrimGalore
```

A textbook benchmark,
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
bowtie2 -q --local -x hg19 -U SRR3177718_trimmed.fq.gz > SRR3177718.bam
samtools sort < SRR3177718.bam > SRR3177718.sort.bam
samtools index SRR3177718.sort.bam
# https://github.com/ncrnalab/agotron_detector
python UCSC_intron_retriever.py | python analyzer.py -g hg19.fa | Rscript annotater.R
```
