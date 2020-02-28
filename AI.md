# Artificial intelligence

## tensorflow

The tensorflow repository is here, https://github.com/tensorflow/tensorflow, and it is relatively easy to install via pip,
```bash
pip install tensorflow
python <<END
import tensorflow as tf
hello = tf.constant('Hello, TensorFlow!')
sess = tf.Session()
print(sess.run(hello))
END
```
Follow https://github.com/aymericdamien/TensorFlow-Examples for readily adaptible examples.

Also

https://github.com/apress/pro-deep-learning-w-tensorflow

## pytorch

The home page is https://pytorch.github.io, and the repository itself https://github.com/pytorch/ with https://github.com/pytorch/examples.

## tensorQTL

AI-derived implementation.

https://github.com/broadinstitute/tensorqtl, https://github.com/broadinstitute/SignatureAnalyzer-GPU

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
# Parquet file
# arrow::install_arrow() to install required runtime libraries
R -e "df <- arrow::read_parquet('GEUVADIS.445_samples.cis_qtl_pairs.chr18.parquet');head(df)"
```
A command-line counterpart is as follows,
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

Taylor-Weiner et al (2019). Scaling computational genomics to millions of individuals with GPUs. *Genome Biol* 20:228,
https://doi.org/10.1186/s13059-019-1836-7
