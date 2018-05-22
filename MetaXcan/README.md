## Issues with more recent version of Python 2.7

The issue was raised to the [MetaXcan GitHub repository](https://github.com/hakyimlab/MetaXcan) for

* Ubuntu 18.04 and Python 2.7.15r1
* Fedora 27 and Python 2.7.15

such that logging.getLogger() was not found.

This was due to confusion between Logging and Python module logging, and fixed by renaming Logging.py to myLogging.py and then adjusting the call from Logging to 
myLogging in M03_betas.py, M04_zscores.py and MetaXcan.py, etc.

It looks the recent version of Python is stricter, which is somewhat expected as with most other compilers. Similar issues were raised while maintaining R packages
for complaints from g++ 8.xx (to be shipped with Fedora 28) which is otherwise OK with g++ 7.x.x.

## Use of the latest databases

While it is possible to use the web interface, https://cloud.hakyimlab.org/user_main, to achieve greater flexibility, the latest databases can be downloaded locally
from [PredictDB Data Repository](http://predictdb.org/).

For instance with [GTEx-V7_HapMap-2017-11-29.tar.gz](https://s3.amazonaws.com/predictdb2/GTEx-V7_HapMap-2017-11-29.tar.gz), we can do the following steps,
```{bash}
mkdir GTEx-V7_HapMap-2017-11-29
cd GTEx-V7_HapMap-2017-11-29
wget https://s3.amazonaws.com/predictdb2/GTEx-V7_HapMap-2017-11-29.tar.gz
tar xvfz GTEx-V7_HapMap-2017-11-29.tar.gz
```
and adjust for the documentation example [test.sh](test.sh) in [V7.sh](V7.sh) as follows,
```{python}
./MetaXcan.py \
--model_db_path /home/jhz22/D/genetics/hakyimlab/ftp/GTEx-V7_HapMap-2017-11-29/gtex_v7_Brain_Amygdala_imputed_europeans_tw_0.5_signif.db \
--covariance /home/jhz22/D/genetics/hakyimlab/ftp/GTEx-V7_HapMap-2017-11-29/gtex_v7_Brain_Amygdala_imputed_eur_covariances.txt.gz \
--gwas_folder data/GWAS \
--gwas_file_pattern ".*gz" \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column BETA \
--pvalue_column P \
--output_file results/V7.csv
```

## Examining weights and related information

[PredictDB FAQs](http://predictdb.org/FAQ.html) point to [a utility in PrediXcan](https://github.com/hakyimlab/PrediXcan/blob/master/Software/query-db.Rmd) for
query, however it is handy to use sqlite3 directory as has been demonstrated in my [TWAS-pipeline](https://github.com/jinghuazhao/TWAS-pipeline). In this case,
we have a utility called [query-db.sql](query-db.sql), which can be called as follows,
```{bash}
sqlite3 gtex_v7_Brain_Amygdala_imputed_europeans_tw_0.5_signif.db < query-db.sql
```
and the weights and extra information are available from files weights.txt and extra.txt, respectively.
