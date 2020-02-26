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
set up virtual environment and install
virtualenv venv
source venv/bin/activate
pip install -r install/requirements.txt .
```

Taylor-Weiner et al (2019). Scaling computational genomics to millions of individuals with GPUs. *Genome Biol* 20:228,
https://doi.org/10.1186/s13059-019-1836-7
