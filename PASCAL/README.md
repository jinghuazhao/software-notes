## Issues with OpenBLAS

When there is issue with xianyi-OpenBLAS-v0.2.12-0-g7e4e195.zip shipped with [PASCAL.zip](http://www2.unil.ch/cbg/images/3/3d/PASCAL.zip), it is recommended to use the [GitHub version](https://github.com/dlampart/Pascal),
as described in [vdi.md](https://github.com/jinghuazhao/GDCT/blob/master/vdi.md) or
```{bash}
git clone https://github.com/xianyi/OpenBLAS
```

## Change of settings.txt

This is necessary since by default pathway analysis is disabled.

Again we use the BMI summary statistics from GIANT,
```{bash}
wget -qO- http://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz | \
gunzip -c | cut -f1,7 | awk -vFS="\t" -vOFS="\t" '(NR>1)' > BMI.pval

Pascal --pval=BMI.pval

```
