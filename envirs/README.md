# Software environments

There are multiple routes to install particular software; one may prefer to install them as standable but it may also come handy use mini-environments such as Anaconda, Miniconda, Linuxbrew or those already in system (e.g. Ubuntu) archive.

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

## Anaconda

Once installed, it is customary to make several channels accessible,

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```
Package in conda-forge include boost, django, glpk, gnuplot, go, gperf, hdf5, ipython, jquery, julia, jupyter, keras, limix, mercurial, miktex, mysql, nano, numpy, pandas, sage, scikit-learn, zlib. Packages in bioconda includes amos, bcftools, beagle, bedops, bedtools, blast, bowtie, bowtie2, bwa, chromhmm, circos, deeptools, emmix, ensembl-vep, fastlmm, fastqc, gatk, gatk4, hclust2, himmer, himmer2, hisat2, igv, impute2, lofreq, mapsplice, mrbayes, ms, nanostat, paml, pbgzip, phylip, picard, plink, plink2, r-wgcna, rsem, rtg-tools, sambamba, samtools, seqkt, sequana, snpeff, snpsift, sra-tools, star, stringtie, tabix, tophat, ucsc-blat, ucsc-liftover, vcftools.

For instance, to install `intervaltree` as required by depict, the following is sufficience,
```bash
conda install intervaltree
```
All the packages installed can be seen with `conda list`. To install java, run following command
```
conda install -c anaconda openjdk
```
Other installations include perl, R. Note that conda under Windows is in typically D:/ProgramData/Anaconda2/Library/bin. Altogether we really need to
```
set path=%path%;D:/ProgramData/Anaconda2;D:/ProgramData/Anaconda2/Library/bin
```
Miniconda is available from https://conda.io/miniconda.html. 

## GitHub

To extract code from GitHub markdown, we do this,

```bash
sudo apt install npm
sudo npm install -g codedown
cat README.md | codedown bash
```
for code in Bash.

**GitKraken** is avaialeble from https://www.gitkraken.com/, e.g., https://release.gitkraken.com/linux/gitkraken-amd64.tar.gz. It can be facilitated with
```bash
sudo apt install libgnome-keyring-common
sudo apt install libgnome-keyring-dev
```

**SmartGit** is available from https://www.syntevo.com/smartgit/, e.g., 
```bash
wget -qO https://www.syntevo.com/downloads/smartgit/smartgit-linux-18_1_4.tar.gz | tar fvxz -
cd smartgit
ln -s $PWD/bin/smartgit.sh $HOME/bin/smartgit.sh
```

**Git-Cola**, https://git-cola.github.io/, can be installed with `sudo apt install git-cola`.

## hg

It is the executable file for Mercurial source code management system,
```bash
sudo apt install mercurial
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

## Parallel computing

It is relevant to have knowledge about GNU parallel and sge. Under Ubuntu, parallel is easily installed as follows,
```{bash}
sudo apt install parallel
```
see also descriptions in other pipelines here. It is perhaps more demanding with sge, e.g., https://peteris.rocks/blog/sun-grid-engine-installation-on-ubuntu-server/.
Example use of slurm can be seen from https://github.com/statgen/SLURM-examples.

## Ubuntu archive

It grows over time, see http://archive.ubuntu.com/ubuntu/pool/universe, including beagle, eigensoft, plink, plink-1.9, among others, which can be installed canonically with ```sudo apt install```.

## Visual Studio Code

There is a pointer from https://github.com/Microsoft/vscode to https://code.visualstudio.com/Download. Once downloaded, it can be installed with
```bash
sudo dpkg -i code_1.23.1-1525968403_amd64.deb
```
but it requires `libgconf-2-4`; when failed to install use `sudo apt --fix-broken install`.

---

## Java

The IDE of choice is NetBeans (e.g., DEPICT and JAM); however 8.1 from `apt install` under Ubuntu 18.04 crashes
so it is suggested to download directly from https://netbeans.org/downloads/. To enable JDK it is helpful to specify `--javahome` option.
```bash
sudo ./netbeans-8.2-linux.sh --javahome /usr/lib/jvm/java-8-oracale
```
or start with `netbeans --javahome /usr/lib/jvm/java-8-oracle` (more convenient to set `alias netbeans='netbeans --javahome /usr/lib/jvm/java-8-oracle'` at `.bashrc`).

NetBeans 9.0 is currently available from https://netbeans.apache.org/download/nb90/; the .zip file can be downloaded and unpacked for use.

For software such as `cutadapt` cython is required,
```bash
sudo apt install cython
```

## Perl
```bash
sudo perl -MCPAN -e shell
install DBI
```
for instance, as used in [VEP](../VEP).

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

PyStan is available with `pip install pystan` which uses matplotlib, https://github.com/matplotlib and Tkinter, established with `sudo apt install python-tk` or `sudo apt install python3-tk`.

```python
import pystan

schools_code = """
data {
    int<lower=0> J; // number of schools
    real y[J]; // estimated treatment effects
    real<lower=0> sigma[J]; // s.e. of effect estimates
}
parameters {
    real mu;
    real<lower=0> tau;
    real eta[J];
}
transformed parameters {
    real theta[J];
    for (j in 1:J)
    theta[j] = mu + tau * eta[j];
}
model {
    eta ~ normal(0, 1);
    y ~ normal(theta, sigma);
}
"""

schools_dat = {'J': 8,
               'y': [28,  8, -3,  7, -1,  1, 18, 12],
               'sigma': [15, 10, 16, 11,  9, 11, 10, 18]}

sm = pystan.StanModel(model_code=schools_code)
fit = sm.sampling(data=schools_dat, iter=1000, chains=4)

import matplotlib.pyplot as plt
def plotGraph():
      fig = fit.plot()
      # plt.show()
      # use the save button or the following command,
      # f.savefig("foo.pdf", bbox_inches='tight')
      return fig

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('foo.pdf')
f = plotGraph()
pp.savefig(f)
pp.close()
```

## Miscellaneous notes
```bash
sudo apt-get install libcanberra-gtk3-module
```

## R

Information on R and RStudio can be seen here, https://github.com/jinghuazhao/Computational-Statistics.

### install.packages and install_github

```r
install.packages("ggplot2",INSTALL_opts="--library=/usr/local/lib/R/site-library/")
install_github("MRCIEU/TwoSampleMR",args="--library=/usr/local/lib/R/site-library",force=TRUE)
```
both supposedly install package to the dedicated location; however this is not always the case and an alternative is to use
```bash
sudo R CMD INSTALL <package_version.tar.gz> -l $R_LIBS
```
to install <package_version.tar.gz> into $R_LIBS.

### The old `tidy.R`

It works as follows,
```bash
function tidy()
{
  export input=$1
R --vanilla <<END
  options(keep.source = FALSE)
  input <- Sys.getenv("input")
  source(input)
  dump(ls(all = TRUE), file = paste0(input,"_out"))
END
}
tidy myfile.R
```
