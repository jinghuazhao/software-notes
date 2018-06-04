# R

## Windows

It seems the --arch x84 option is very useful for using all available RAM; to make sure use call such as `D:\Program Files\R\R-3.5.0\bin\x64\R.exe"`.

When this fails, remove large objects in your code and start R with `--vanilla` option.

## Fedora 28

The guest additions is furnished with
```bash
sudo dnf install gcc kernel-devel kernel-headers dkms make bzip2 perl
cd /run/media/jhz22/VBox_GAs_5.2.12/
sudo ./VBoxLinuxAdditions.run
```

## Ubuntu 18.04

### R installation

```{bash}
sudo apt install r-base-core
sudo apt install r-base-dev
```
and R_LIBS is set from .bashrc
```{bash}
export R_LIBS=/usr/local/lib/R/site-library/
```
Note that in fact `html.start()` in R points to /usr/local/lib/R/library/ instead, see below example in `MendelianRandomization`.

### R-devel

Under Fedora 28, the following are necessary,
```bash
sudo dnf install gcc-c++
sudo dnf install gcc-gfortran
sudo dnf install compat-gcc-34-g77
sudo dnf install java-openjdk-devel
sudo dnf install pcre-devel
sudo dnf install readline-devel
sudo dnf install libcurl-devel
sudo dnf install libX11 libX11-devel libXt libXt-devel
sudo dnf install bzip2-devel
sudo dnf install xz-devel
sudo dnf install texlive-collection-latex
sudo dnf install texlive-collection-fontsextra
sudo dnf install texinfo-tex
./configure
```

### RStudio

The distribution has problem loading or creating R script, so it is tempting to install from https://github.com/rstudio/rstudio/. This involves running scripts under directory dependencies/, 
```{bash}
./install-dependencies-debian --exclude-qt-sdk
```
and then the following steps,
```{bash}
mkdir build
cd build
cmake .. -DRSTUDIO_TARGET=Desktop -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local/lib/rstudio
```
However, there is error with Java and Java 8 is required, see https://tecadmin.net/install-oracle-java-8-ubuntu-via-ppa/.
```{bash}
sudo add-apt-repository ppa:webupd8team/java
sudo apt-get update
sudo apt-get install oracle-java8-installer
sudo apt-get install oracle-java8-set-default
java -version
```
However, compile error is still persistent except when dropping the option `--exclude-qt-sdk` but unloadable.

It is therefore recommended to get around with RStudio daily builds, https://dailies.rstudio.com/.
