# ldsc

## Partitioned heritability

The [wiki documentation](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability) script really should be as follows,
```bash

export GIANT_BMI=GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt

setup()
{
  if [ ! -f $GIANT_BMI ]; then
     wget http://portals.broadinstitute.org/collaboration/giant/images/b/b7/$GIANT_BMI.gz
     gunzip $GIANT_BMI.gz
  fi
  if [ ! -f w_hm3.snplist ]; then
     wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
     bzip2 -d w_hm3.snplist.bz2
  fi
  python munge_sumstats.py --sumstats $GIANT_BMI --merge-alleles w_hm3.snplist --out BMI --a1-inc
}

# It fails to munge, so we use brute-force. For P-to-Z implementation in C/C++, see
# https://stackoverflow.com/questions/27830995/inverse-cumulative-distribution-function-in-c
# https://stackoverflow.com/questions/22834998/what-reference-should-i-use-to-use-erf-erfc-function

awk 'NR>1' $GIANT_BMI > 1
awk 'NR>1' w_hm3.snplist | sort -k1,1 | join -j1 1 - | awk -f CLEAN_ZSCORES.awk > BMI.sumstats
R --vanilla -q <<END
BMI <- read.table("BMI.sumstats",col.names=c("SNP","A1","A2","Z","N"))
BMI <- within(BMI, {Z=sign(Z)*qnorm(abs(Z)/2)})
z <- gzfile("BMI.sumstats.gz","w")
write.table(BMI,file=z,quote=FALSE,row.names=FALSE)
close(z)
END
```
where we use [CLEAN_ZSCORES.awk](CLEAN_ZSCORES.awk) to align SNPs between sumstats and reference.

Now the partition heritability and cell-type group analysis proceed as follows,
```bash
python ldsc.py --h2 BMI.sumstats.gz\
        --ref-ld-chr baseline_v1.1/baseline.\
        --w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.\
        --overlap-annot\
        --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC.\
        --out BMI_baseline

python ldsc.py --h2 BMI.sumstats.gz\
        --w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.\
        --ref-ld-chr 1000G_Phase3_cell_type_groups/cell_type_group.3.,baseline_v1.1/baseline.\
        --overlap-annot\
        --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC.\
        --out BMI_CNS\
        --print-coefficients
```
NB it is assumed that [all the required data](https://data.broadinstitute.org/alkesgroup/LDSCORE/) have been made available.

## Test

The mysterious commands shown in the wiki documentation are actually realised after this,
```
sudo apt install python-nose
```
