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
ln -s /genetics/bin/.vep
```
for instance.
