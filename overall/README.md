# Overall setup

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

rem Windows
path D:\Program Files\Oracle\VirtualBox
VBoxManage modifyhd --compact "ubuntu18.04.vdi"
```

[vdi.md](https://github.com/jinghuazhao/GDCT/blob/master/vdi.md) as in GWAS-2017 and now listed in [GDCT](https://github.com/jinghuazhao/GDCT)

Since one may allocate only part of RAM to VirtualBox, it is often necessary to run program under MS-DOS, e.g., sections on DEPICT.

## Visual Studio Code

There is a pointer from https://github.com/Microsoft/vscode to https://code.visualstudio.com/Download. Once downloaded, it can be installed with
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

## [R and RStudio](../R)

## Perl
```bash
sudo perl -MCPAN -e shell
install DBI
```
for instance, as used in [VEP](../VEP).

## Java

The IDE of choice is NetBeans (e.g., DEPICT and JAM); however 8.1 from `apt install` under Ubuntu 18.04 crashes
so it is suggested to download directly from https://netbeans.org/downloads/. To enable JDK it is helpful to specify `--javahome` option.
```bash
sudo ./netbeans-8.2-linux.sh --javahome /usr/lib/jvm/java-8-oracale
```
or start with `netbeans --javahome /usr/lib/jvm/java-8-oracle` (more convenient to set `alias netbeans='netbeans --javahome /usr/lib/jvm/java-8-oracle'` at `.bashrc`).

For software such as `cutadapt` cython is required,
```bash
sudo apt install cython
```
## Python

To install a particular version of package, e.g.,
```bash
sudo -H pip install pandas==0.20.1
```
which is required by DEPICT's `munge_sumstats.py`. Other pip options include `uninstall`.

The python programs in [agotron_detector](https://github.com/ncrnalab/agotron_detector) requires MySQL and can be installed as follows,
```bash
sudo apt-get install python-dev libmysqlclient-dev
sudo pip install MySQL-python
```

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

## JAGS-4.3.0

These are required at least,
```bash
sudo dnf install automake
sudo dnf install lapack-devel
sudo dnf install mercurial
```

## Linuxbrew

Follow http://linuxbrew.sh/ and possibly https://docs.brew.sh
```bash
sudo apt-get install build-essential
sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
echo 'export PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"' >>~/.profile
echo 'export MANPATH="/home/linuxbrew/.linuxbrew/share/man:$MANPATH"' >>~/.profile
echo 'export INFOPATH="/home/linuxbrew/.linuxbrew/share/info:$INFOPATH"' >>~/.profile
PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"
```

## hg

It is the executable file for Mercurial source code management system,
```bash
sudo apt install mercurial
```

## Ubuntu lists

It grows over time, see http://archive.ubuntu.com/ubuntu/pool/universe, including beagle, plink, plink-1.9, among others.


## Miscellaneous notes
```bash
sudo apt-get install libcanberra-gtk3-module
```
