# LDpred

Now it is rather simple to install, e.g.,
```bash
pip install ldpred
```
or
```bash
pip install --user ldpred
```
Nevertheless it is rather laborous to read through the documentation and we try the GitHub version as well.
```bash
git clone https://github.com/bvilhjal/ldpred
cd ldpred
cd ldpred
python test.py
```
rather than a Bash version from scratch which is also possible,
```bash
cd ldpred
export test=ZDNCEO
python coord_genotypes.py --gf=../test_data/LDpred_data_p0.001_train_0 \
                          --vgf=../test_data/LDpred_data_p0.001_test_0 \
                          --ssf=../test_data/LDpred_data_p0.001_ss_0.txt \
                          --N=10000 \
                          --out=$test.coord.hdf5
python LDpred.py --coord=$test.coord.hdf5  \
                 --ld_radius=100 \
                 --local_ld_file_prefix=$test \
                 --PS=0.001 \
                 --N=10000  \
                 --out=$test
python validate.py --vgf=../test_data/LDpred_data_p0.001_test_0 \
                   --rf=$test \
                   --out=$test
```
