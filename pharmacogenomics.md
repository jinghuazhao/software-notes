
### R/Pi

```r
devtools::install_bioc("Pi",build_vignettes = TRUE)
library(Pi)
browseVignette("Pi")
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

Reference

Fang H; ULTRA-DD Consortium, De Wolf H, Knezevic B, Burnham KL, Osgood J, Sanniti A, Lledó Lara A, Kasela S, De Cesco S, Wegner JK, Handunnetthi L, McCann FE, Chen L, Sekine T, Brennan PE, Marsden BD, Damerell D, O'Callaghan CA, Bountra C, Bowness P, Sundström Y, Milani L, Berg L, Göhlmann HW, Peeters PJ, Fairfax BP, Sundström M, Knight JC.
A genetics-led approach defines the drug target landscape of 30 immune-related traits.
at Genet. 2019 Jul;51(7):1082-1091. doi: 10.1038/s41588-019-0456-1. Epub 2019 Jun 28.
