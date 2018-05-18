# PyLMM

The software is rare with its setup for GEI studies accounting for polygenic effects.

It is necessary to do the following steps to get it going,

1. git clone https://github.com/nickFurlotte/pylmm
2. change
  * build/scripts-2.7/pylmmGWAS.py, line 207, from `keep = True - v` to `keep = True ^ v`
  * pylmm/lmm.py, line 189, from `if X0 == None:` to `if X0.all == None:`; line 193, from `x = True - np.isnan(Y)` to `x = True ^ np.isnan(Y)`; line 272, from `if X == None: X = self.X0t` to `if X.all == None: X = self.X0t`.
  * build/scripts-2.7/input.py, line 190, from `x = True ^ np.isnan(G)` to `x = True ^ np.isnan(G)`.
3. python setup.py install
4. bash [test.sh](test.sh) which is the documentation call,
```{bash}
pylmmGWAS.py -v --bfile data/snps.132k.clean.noX --kfile data/snps.132k.clean.noX.pylmm.kin --phenofile data/snps.132k.clean.noX.fake.phenos out.foo
```
