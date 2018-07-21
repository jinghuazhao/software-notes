# Miscellaneous software

The Ubuntu archive, http://archive.ubuntu.com/ubuntu/pool/universe, is installed canonically with ```sudo apt install```. The source code is also available from the archive. We can trick the test even with compiling errors, e.g., 
```bash
sudo apt install loki
wget -qO- http://archive.ubuntu.com/ubuntu/pool/universe/l/loki/loki_2.4.7.4.orig.tar.gz | tar fvxz -
cd loki
./configure
# there are errors in compilng
make
cp /usr/bin/prep presrc
cp /usr/bin/loki lokisrc
cd test
# we are actually fine
make
```
There is a variety of packages in bioconda, https://bioconda.github.io/recipes.html# and conda-forge, https://conda-forge.org/feedstocks/.

## allegro

This is according to https://www.decode.com/software/,
```bash
sudo `which conda` install allegro
```
> Gudbjartsson DF, Thorvaldsson T, Kong A, Gunnarsson G, Ingolfsdottir A (2005). Allegro version 2. *Nature Genetics* 37:1015–1016

https://www.nature.com/articles/ng1005-1015?foxtrotcallback=true

## beast2-mcmc, beast2-mcmc-doc

It is Bayesian MCMC phylogenetic inference.
```bash
sudo apt install beast2-mcmc
```

## fastlink, fastlink-doc

It is the old merry fastlink 4.1P.
```bash
sudo apt install fastlink
```

## loki, loki-doc

It implements MCMC method for linkage analysis.
```bash
sudo apt install loki
```
See usr/share/doc/loki-doc/ for documentation.

## merlin

```bash
sudo `which conda` install merlin
```

## simwalk2

```bash
sudo `which conda` install simwalk2
```
