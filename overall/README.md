# Overall setup

Contents in this section are generic, therefore worthwhile to note and maintain for a variety of applications.

## Resources

https://www.genomicsengland.co.uk/ and [centres](https://www.genomicsengland.co.uk/taking-part/genomic-medicine-centres/)

http://www.ukbiobank.ac.uk/

https://www.astrazeneca.com/

https://opensource.microsoft.com/

## Imputation

Sanger Imputation Service, https://imputation.sanger.ac.uk/

The Michigan imputation server, https://imputationserver.sph.umich.edu/start.html

Code for 1000Genomes populations, http://www.internationalgenome.org/category/population/

## Oracle VirtualBox

To set up shared folders and enforce shared clipboard for bidirectional copy between Linux and Windows,
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

[vdi.md](https://github.com/jinghuazhao/GDCT/blob/master/vdi.md) as in GWAS-2017 and now listed in [GDCT](https://github.com/jinghuazhao/GDCT)

Since one may allocate only part of RAM to VirtualBox, it is often necessary to run program under MS-DOS, e.g., sections on DEPICT.

## Visual Studio Code

There is a pointer from https://github.com/Microsoft/vscode to https://code.visualstudio.com/Download. Once download, it can be installed with
```bash
sudo dpkg -i code_1.23.1-1525968403_amd64.deb
```
but it requires `libgconf-2-4`; when failed to install use `sudo apt --fix-broken install`.

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

## Java

The IDE of choice is NetBeans (e.g., DEPICT and JAM); however 8.1 from `apt install` under Ubuntu 18.04 crashes
so it is suggested to download directly from https://netbeans.org/downloads/.

## Anaconda

To install intervaltree as required by depict,
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```
followed by `conda install intervaltree`.

All the available packages can be seen with `conda list`. To install java, run following command
```
conda install -c anaconda openjdk
```
Note that conda under Windows is in typically D:/ProgramData/Anaconda2/Library/bin. Altogether we really need to
```
set path=%path%;D:/ProgramData/Anaconda2;D:/ProgramData/Anaconda2/Library/bin
```
