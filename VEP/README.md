# VEP

The description is available from http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html.
```bash
git clone https://github.com/Ensembl/ensembl-vep.git

cd ensembl-vep
git pull
git checkout release/92
perl INSTALL.pl
```
The last line requires modules DBI, Build as described in [Overall](../overall).

Lastly, VEP requires .vep directory at $HOME which can be derived from a centrally-installed VEP under Linux,
```bash
cd $HOME
ln -s /genetics/ensembl-vep/.vep
```
for instance.

It is slow to get those databases, so one may prefer to get them directly from ftp://ftp.ensembl.org/pub/release-92/variation/VEP/ and unpack into the .vep directory.
