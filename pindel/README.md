# pindel

The software can be obtained from https://github.com/genome/pindel.

After htslib is installed, the canonical instruction is to issue
```bash
git clone https://github.com/samtools/htslib
cd htslib
make
sudo make install
cd -
git clone https://github.com/genome/pindel
cd pindel
./INSTALL ../htslib
```
It is `standard' to have complaints about pindel.cpp, bddate.cpp and genotyping.cpp,
for `abs()` rather than `fabs()` from the header file `cmath` have been used. The
issue goes away when `abs` is replaced with `fabs` and in the case of bddata.cpp,
it is also necessary to invoke the header, i.e.,
```cpp
#include <cmath>
```
