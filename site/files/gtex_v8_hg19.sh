#!/usr/bin/bash

# https://github.com/hmgu-itg/VCF-liftover
# rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz .

for f in gtex_v8.eqtl_annot_rsid gtex_v8.eqtl_annot
do
  export f=${f}
  tabix -f ${f}.vcf.gz
  (
    seq 22 | \
    parallel -j1 --env f -C' ' '
      vcf-liftover.sh hg38ToHg19.over.chain ${f}.vcf.gz ${f}.hg19-chr{}.vcf.gz chr{} dontstop fast
      zcat ${f}.hg19-chr{}.vcf.gz
      rm ${f}.hg19-chr{}.vcf.gz
    '
  ) | \
  bgzip -f > ${f}.hg19.vcf.gz
done

(
  seq 22 | \
  parallel -j1 -C' ' 'cat hg38ToHg19.over.chain.chr{}.offset'
) > hg38ToHg19.over.chain.offset
rm hg38ToHg19.over.chain.chr*.offset

chr=chr1
pos=1234567
file=hg38ToHg19.over.chain.offset
offset=$(awk -v c=$chr -v p=$pos '$1==c && p>=$2 && p<=$3{print $4}' $file)
echo $(( $pos + $offset ))

zgrep -v -e Written -e Lifting gtex_v8.eqtl_annot_rsid.hg19.vcf.gz | \
bgzip -f > gtex_v8.eqtl.hg19.vcf.gz

gunzip -c gtex_v8.eqtl.hg19.vcf.gz | \
cut -f6 | tr '||' '\n' | tr '@' ' ' | tr '=' ' '| cut -d' ' -f2 | sort | uniq > Tissues
