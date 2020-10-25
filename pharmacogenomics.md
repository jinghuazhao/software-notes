### Python/chembl_webresource_client

Web: https://github.com/chembl/chembl_webresource_client
```bash
pip install chembl_webresource_client
```
UniProtID mapping API: https://www.uniprot.org/help/api_idmapping
```python3
import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC+ID',
'to': 'CHEMBL_ID',
'format': 'tab',
'query': 'P40925 P40926 O43175 Q9UM73 P97793'
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
print(response.decode('utf-8'))
``` 

### R/Pi

```r
devtools::install_bioc("Pi",build_vignettes = TRUE)
library(Pi)
browseVignettes("Pi")
```
Note that `S4Vectors` may conflict with the R-devel configurations.

The figshare database as described in the paper can be examined as follows,

```bash
gunzip pi_database.sql.gz
mysql -p -u $USER -e "create database pi;"
mysql -p -u $USER pi < pi_database.sql
mysql -p -u $USER pi <<END
show tables;
desc pi_priority;
desc pi_trait;
desc pi_genomic;
desc pi_drug;
desc pi_category;
desc pi_domain;
desc pi_pdb;
select * from pi_priority where trait='RA' and rank<=150 order by rank;
select * from pi_priority where trait='RA' and crosstalk_node='Y' order by rank;
END
```
It is relatively easier to convert these as an [R object](files/pi_database.rda),
```r
library(RMySQL)
user <- "my user name"
password <- "mypassword"
mydb = dbConnect(MySQL(), user=user, password=password, dbname='pi')
tbllist <- dbListTables(mydb)
pi_category <- dbGetQuery(mydb, paste0("select * from pi_category;"))
pi_domain <- dbGetQuery(mydb, paste0("select * from pi_domain;"))
pi_drug <- dbGetQuery(mydb, paste0("select * from pi_drug;"))
pi_genomic <- dbGetQuery(mydb, paste0("select * from pi_genomic;"))
pi_pdb <- dbGetQuery(mydb, paste0("select * from pi_pdb;"))
pi_priority <- dbGetQuery(mydb, paste0("select * from pi_priority;"))
pi_trait <- dbGetQuery(mydb, paste0("select * from pi_trait;"))
# target category
head(pi_category)
# target superfamily druggable
head(pi_domain)
# target uniprot pdb_chain pdb chain pocket
head(pi_pdb)
# trait target rank rating nGene cGene eGene seed fGene pGene dGene gwas crosstalk_node num_neighbor approved phased druggable_category druggable_domain num_pdb num_pdb_with_druggable_pocket magnitude direction description
head(pi_priority)
traits <- c("AS","CRO","IGE","MS","RA","T1D","UC")
# trait target max_phase drug mechanism_of_action action_type source
subset(pi_drug,trait%in%traits)
# trait target type name snp snp_type pvalue
subset(pi_genomic,trait%in%traits)
save(pi_category,pi_domain,pi_drug,pi_genomic,pi_pdb,pi_priority,pi_trait,file="pi_database.rda")
```

### R/dbparser

https://cran.r-project.org/web/packages/dbparser/

### R/rDGIdb

https://www.bioconductor.org/packages/release/bioc/html/rDGIdb.html

http://www.dgidb.org/


### R/sunburstR

https://github.com/timelyportfolio/sunburstR

Techinical aspect: https://stackoverflow.com/questions/12926779/how-to-make-a-sunburst-plot-in-r-or-python

Example: https://bl.ocks.org/kerryrodden/7090426

## Reference

Cotto KC, Wagner AH, Feng YY, Kiwala S, Coffman AC, Spies G, Wollam A, Spies NC, Griffith OL, Griffith M. DGIdb 3.0: a redesign and expansion of the drug-gene interaction database. *Nucleic Acids Res.* 2018 Jan 4;46(D1):D1068-D1073. doi: 10.1093/nar/gkx1143. PMID: 29156001; PMCID: PMC5888642.

Fang H; ULTRA-DD Consortium, De Wolf H, Knezevic B, Burnham KL, Osgood J, Sanniti A, Lledó Lara A, Kasela S, De Cesco S, Wegner JK, Handunnetthi L, McCann FE, Chen L, Sekine T, Brennan PE, Marsden BD, Damerell D, O'Callaghan CA, Bountra C, Bowness P, Sundström Y, Milani L, Berg L, Göhlmann HW, Peeters PJ, Fairfax BP, Sundström M, Knight JC.
A genetics-led approach defines the drug target landscape of 30 immune-related traits.
*Nat Genet*. 2019 Jul;51(7):1082-1091. doi: 10.1038/s41588-019-0456-1. Epub 2019 Jun 28.
