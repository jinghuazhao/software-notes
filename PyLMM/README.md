# PyLMM

The software is rare with its setup for GEI studies accounting for polygenic effects.


## pylmm

It is necessary to do the following steps to get it going,

1. git clone https://github.com/nickFurlotte/pylmm
2. make the following changes,
  * build/scripts-2.7/pylmmGWAS.py, line 207, from `keep = True - v` to `keep = True ^ v`
  * pylmm/lmm.py, line 189, from `if X0 == None:` to `if X0.all == None:`; line 193, from `x = True - np.isnan(Y)` to `x = True ^ np.isnan(Y)`; line 272, from `if X == None: X = self.X0t` to `if X.all == None: X = self.X0t`.
  * build/scripts-2.7/input.py, line 190, from `x = True ^ np.isnan(G)` to `x = True ^ np.isnan(G)`.
3. python setup.py install
4. bash [test.sh](test.sh) which is the documentation call,
```{bash}
pylmmGWAS.py -v --bfile data/snps.132k.clean.noX --kfile data/snps.132k.clean.noX.pylmm.kin --phenofile data/snps.132k.clean.noX.fake.phenos out.foo
```

## pylmm_zarlab

Make the following changes to lmm and lmmGWAS similar to pylmm, and then issue `bash run_tests.sh`.

```
EReading SNP input...
Read 1219 individuals from data/snps.132k.clean.noX.fam
Reading kinship...
Read the 1219 x 1219 kinship matrix in 1.139s 
1 number of phenotypes read
Traceback (most recent call last):
  File "scripts/pylmmGWAS.py", line 308, in <module>
    keep = True - v
TypeError: numpy boolean subtract, the `-` operator, is deprecated, use the bitwise_xor, the `^` operator, or the logical_xor function instead.
EReading PLINK input...
Read 1219 individuals from data/snps.132k.clean.noX.fam
Traceback (most recent call last):
  File "scripts/pylmmKinship.py", line 127, in <module>
    K_G = lmm.calculateKinshipIncremental(IN, numSNPs=options.numSNPs,
AttributeError: 'module' object has no attribute 'calculateKinshipIncremental'
EReading PLINK input...
Read 1219 individuals from data/snps.132k.clean.noX.fam
Traceback (most recent call last):
  File "scripts/pylmmKinship.py", line 127, in <module>
    K_G = lmm.calculateKinshipIncremental(IN, numSNPs=options.numSNPs,
AttributeError: 'module' object has no attribute 'calculateKinshipIncremental'
E
======================================================================
ERROR: test_GWAS (tests.test_lmm.test_lmm)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/tests/test_lmm.py", line 41, in test_GWAS
    TS,PS = lmm.GWAS(Y,snps,K,REML=True,refit=True)
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/pylmm/lmm.py", line 192, in GWAS
    L = LMM(Y, K, Kva, Kve, X0)
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/pylmm/lmm.py", line 301, in __init__
    x = True - np.isnan(Y)
TypeError: numpy boolean subtract, the `-` operator, is deprecated, use the bitwise_xor, the `^` operator, or the logical_xor function instead.

======================================================================
ERROR: test_calculateKinship (tests.test_lmm.test_lmm)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/tests/test_lmm.py", line 25, in test_calculateKinship
    K = lmm.calculateKinship(snps)
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/pylmm/lmm.py", line 135, in calculateKinship
    mn = W[True - np.isnan(W[:, i]), i].mean()
TypeError: numpy boolean subtract, the `-` operator, is deprecated, use the bitwise_xor, the `^` operator, or the logical_xor function instead.

======================================================================
ERROR: test_pylmmGWASScript (tests.test_pylmmGWAS.test_pylmmGWAS)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/tests/test_pylmmGWAS.py", line 24, in test_pylmmGWASScript
    with (open(self._outputFile, 'r')) as ansFile:
IOError: [Errno 2] No such file or directory: 'data/pylmmGWASTestOutput'

======================================================================
ERROR: test_pylmmKinshipScript1 (tests.test_pylmmKinship.test_pylmmKinship)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/tests/test_pylmmKinship.py", line 24, in test_pylmmKinshipScript1
    K = np.fromfile(open(self._outputFile, 'r'), sep=" ")
IOError: [Errno 2] No such file or directory: 'data/pylmmKinshipTestOutput'

======================================================================
ERROR: test_pylmmKinshipScript2 (tests.test_pylmmKinship.test_pylmmKinship)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/jhz22/D/genetics/ucla/pylmm_zarlab/tests/test_pylmmKinship.py", line 38, in test_pylmmKinshipScript2
    K = np.fromfile(open(self._outputFile, 'r'), sep=" ")
IOError: [Errno 2] No such file or directory: 'data/pylmmKinshipTestOutput'

----------------------------------------------------------------------
Ran 5 tests in 2.776s

FAILED (errors=5)
```
We can have a test of GxE analysis as this,
```{bash}
sudo python setup.py install
cd pylmm
python pylmm_GXE.py
```
In general, we can see options for GxE analysis from command `pylmmGWAS.py` under bash.

