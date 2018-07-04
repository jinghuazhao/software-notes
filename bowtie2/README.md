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
