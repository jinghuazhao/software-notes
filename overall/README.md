# Overall setup

## Resources

https://www.genomicsengland.co.uk/ and [centres](https://www.genomicsengland.co.uk/taking-part/genomic-medicine-centres/)

http://www.ukbiobank.ac.uk/

https://www.astrazeneca.com/

## Imputation

Sanger Imputation Service, https://imputation.sanger.ac.uk/

The Michigan imputation server, https://imputationserver.sph.umich.edu/start.html

Code for 1000Genomes populations, http://www.internationalgenome.org/category/population/

## Parallel computing

It is relevant to have knowledge about GNU parallel and sge. Under Ubuntu, parallel is easily installed as follows,
```{bash}
sudo apt install parallel
```
see also descriptions in other pipelines here. It is perhaps more demanding with sge, e.g., https://peteris.rocks/blog/sun-grid-engine-installation-on-ubuntu-server/.
Example use of slurm can be seen from https://github.com/statgen/SLURM-examples.

## [R and RStudio](../R).

## Perl
```bash
sudo perl -MCPAN -e shell
install DBI
```
for instance, as used in [VEP](../VEP).

## Anaconda

To install intervaltree as required by depict,
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```
followed by `conda install intervaltree`.

All the available packages can be seen with `conda list`. Note that conda under WIndows is in typically D:/ProgramData/Anaconda2/Library/bin.

## [vdi.md](https://github.com/jinghuazhao/GDCT/blob/master/vdi.md) as in GWAS-2017 and now listed in [GDCT](https://github.com/jinghuazhao/GDCT)

## Oracle VirtualBox

To enforce sharedclipboard for bidirectional copy between Linux and Windows,
```bash
# shared folders
sudo mount -t vboxsf -o uid=jhz22 C /home/jhz22/C
sudo mount -t vboxsf -o uid=jhz22 D /home/jhz22/D
# shared clipboard
killall VBoxClient
sudo VBoxClient-all
```
Here are the steps, quoting http://www.netreliant.com/news/8/17/Compacting-VirtualBox-Disk-Images-Linux-Guests.html, for compressing large .vdi:
```bash
# Linux
dd if=/dev/zero of=zerofillfile bs=1M

# Windows
path D:\Program Files\Oracle\VirtualBox
VBoxManage modifyhd --compact "ubuntu18.04.vdi"
```
