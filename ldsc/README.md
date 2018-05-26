# ldsc

## Partitioned heritability

The [wiki documentation](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability) script really should be as follows,
```bash
# https://github.com/bulik/ldsc
# only keep HapMap3 SNPs

export GIANT_BMI=GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt

setup()
{
  if [ -f w_hm3.snplist ]; then
     wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
     bzip2 -d w_hm3.snplist.bz2
  fi
  python munge_sumstats.py --sumstats $GIANT_BMI --merge-alleles w_hm3.snplist --out BMI --a1-inc
}

# somehow munge_sumstats.py fails under my Windwos 10 and Ubuntu 18.04, so I trick going with P as Z
awk 'NR>1' $GIANT_BMI > 1
awk 'BEGIN{print "SNP","A1","A2","Z","N"}' > BMI.sumstats
awk 'NR>1' w_hm3.snplist | sort -k1,1 | join -j1 1 - | awk -f CLEAN_ZSCORES.awk >> BMI.sumstats
gzip -f BMI.sumstats
rm 1
```
where we use [CLEAN_ZSCORES.awk][CLEAN_ZSCORES.awk] to align SNPs between sumstats and reference.

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
        --ref-ld-chr 1000G_Phase3_cell_type_groups/cell_type_group.3.\
        --overlap-annot\
        --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC.\
        --out BMI_CNS\
        --print-coefficients
```
NB it is assumed that [all the required data](https://data.broadinstitute.org/alkesgroup/LDSCORE/) have been available.

## Test

The mysterious commands shown in the wiki documentation are actually realised after this,
```
sudo apt install python-nose
```
