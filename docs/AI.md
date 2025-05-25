# Artificial intelligence

## DrugAssist

Ye, Geyan, et al. "DrugAssist: A Large Language Model for Molecule Optimization." Briefings in Bioinformatics, vol. 26, no. 1, 2025, p. bbae693. Oxford University Press.

<https://github.com/blazerye/DrugAssist> & <https://huggingface.co/blazerye/DrugAssist-7B>

## nucleotide-transformer

<https://github.com/instadeepai/nucleotide-transformer>

## Pyro

<https://pyro.ai/examples/index.html>

## pytorch

The home page is <https://pytorch.github.io>, and the repository itself <https://github.com/pytorch/> with <https://github.com/pytorch/examples>.

## tensorflow

### 1.x

The tensorflow repository is here, <https://github.com/tensorflow/tensorflow>, and it is relatively easy to install via pip,
```bash
pip install tensorflow
python <<END
import tensorflow as tf
hello = tf.constant('Hello, TensorFlow!')
sess = tf.Session()
print(sess.run(hello))
END
```
Follow <https://github.com/aymericdamien/TensorFlow-Examples> for readily adaptible examples.

Also

<https://github.com/apress/pro-deep-learning-w-tensorflow>

### 2.x

python <<END
import tensorflow as tf

hello = tf.constant('Hello, TensorFlow!')
print(hello.numpy().decode('utf-8'))
END

## tensorQTL

AI-derived implementation.

<https://github.com/broadinstitute/tensorqtl>, <https://github.com/broadinstitute/SignatureAnalyzer-GPU>

```bash
module load python/3.6
git clone git@github.com:broadinstitute/tensorqtl.git
cd tensorqtl
# set up virtual environment and install
virtualenv venv
source venv/bin/activate
pip install -r install/requirements.txt .
cd example
wget https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bed
wget https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bim
wget https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.fam   
wget https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.covariates.txt
wget https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.expression.bed.gz
# Jupyter notebook
sed -i 's/filtered/filtered.nodup/g' tensorqtl_examples.ipynb
# csd3
hostname
jupyter notebook --ip=127.0.0.1 --no-browser --port 8081
# local host
ssh -4 -L 8081:127.0.0.1:8081 -fN hostname.hpc.cam.ac.uk
firefox <generated URL from jupyter notebook command above> &
```
Note that a Parquet file is generated we use SparkR,
```bash
module load spark/2.4.0-bin-hadoop2.7
```
followed by
```r
library(SparkR)
sparkR.session()
df <- read.parquet("GEUVADIS.445_samples.cis_qtl_pairs.chr18.parquet")
head(df)
```
to get
```
> dim(df)
[1] 2927819       9
> head(df)
       phenotype_id          variant_id tss_distance        maf ma_samples
1 ENSG00000263006.6 chr18_10644_C_G_b38       -98421 0.01685393         15
2 ENSG00000263006.6 chr18_10847_C_A_b38       -98218 0.01910112         17
3 ENSG00000263006.6 chr18_11275_G_A_b38       -97790 0.02471910         22
4 ENSG00000263006.6 chr18_11358_G_A_b38       -97707 0.02471910         22
5 ENSG00000263006.6 chr18_11445_G_A_b38       -97620 0.02359551         21
6 ENSG00000263006.6 chr18_13859_G_C_b38       -95206 0.02471910         22
  ma_count pval_nominal       slope  slope_se
1       15    0.5808729 -0.11776078 0.2131254
2       17    0.1428839 -0.29872555 0.2035047
3       22    0.7452308  0.05461900 0.1679810
4       22    0.7452308  0.05461900 0.1679810
5       21    0.6032759  0.08937798 0.1718505
6       22    0.7452308  0.05461900 0.1679810
```
An alternative is to tweak the R package `arrow`.

The command-line counterpart is as follows,
```bash
export plink_prefix_path=GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup
export expression_bed=GEUVADIS.445_samples.expression.bed.gz
export covariates_file=GEUVADIS.445_samples.covariates.txt
export prefix=GEUVADIS.445_samples

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode trans
```
Again one can read the Parquet format output.

Taylor-Weiner et al (2019). Scaling computational genomics to millions of individuals with GPUs. *Genome Biol* 20:228,
<https://doi.org/10.1186/s13059-019-1836-7>
