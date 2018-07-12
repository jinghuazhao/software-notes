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
It could be useful with ``brew reinstall```. See
> Raghavachari N, Garcia-Reyero N (eds.) (2018), Gene Expression Analysis-Methods and Protocols, https://www.springer.com/us/book/9781493978335, Chapter 15, Springer.

Nevertheless it may be slower, e.g., tophat, compared to ```sudo apt install```.

# sra-toolkit, samtools, bamtools, bcftools, tophat, vcftools, blast (ncbi-blast+)

These are very straightforward under Ubuntu, e.g.,
```bash
sudo apt install sra-toolkit
prefetch -v SRR3534842
fastq-dump --split-files --gzip SRR3534842
```
the SRR3534842.sra from prefetch is actually at $HOME/ncbi/public/sra which is split into
`SRR3534842_1.fastq.gz`, `SRR3534842_2.fastq.gz` at the current directory.

See https://www.biostars.org/p/111040/.

To run ```top```, see https://ccb.jhu.edu/software/tophat/tutorial.shtml
```bash
wget https://ccb.jhu.edu/software/tophat/downloads/test_data.tar.gz
tar xvfz test_data.tar.gz
cd test_data
tophat -r 20 test_ref reads_1.fq reads_2.fq
```

# cutadapt, TrimGalore

A prerequesite is to install cython.
```bash
git clone https://github.com/marcelm/cutadapt
cd cutadapt
sudo python setup.py install
git clone https://github.com/FelixKrueger/TrimGalore
```

# A textbook benchmark

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
See also
> Martin J, Schackwitz W, Lipzen A (2018). Genomic Sequence Variation Analysis by Resequencing, in de Vries RP, Tsang A, Grigoriev IV (ed) Fungal Genomics-Methods and Protocols, 2e, Chapter 18, 229-239, Springer.

> Morin PA, Foote AD, Hill CM, Simon-Bouhet B, Lang AR, Louis M (2018). SNP Discovery from Single and Multiplex Genome Assemblies of Non-model Organisms, in Steven R. Head et al. (eds.), Next Generation Sequencing: Methods and Protocols, Chapter 9, 113-144, Springer.

> Raghavachari N, Garcia-Reyero N (eds.) (2018), Gene Expression Analysis-Methods and Protocols, Springer.
