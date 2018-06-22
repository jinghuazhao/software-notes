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

## Resources

https://www.genomicsengland.co.uk/ and [centres](https://www.genomicsengland.co.uk/taking-part/genomic-medicine-centres/)

http://human-phenotype-ontology.github.io/about.html

http://www.ukbiobank.ac.uk/

https://www.astrazeneca.com/

https://opensource.microsoft.com/

## Imputation

Sanger Imputation Service, https://imputation.sanger.ac.uk/

The Michigan imputation server, https://imputationserver.sph.umich.edu/start.html

Code for 1000Genomes populations, http://www.internationalgenome.org/category/population/

