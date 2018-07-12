# 12-7-2018 JHZ

git add README.md
git commit -m "README"
git add topics.md
git commit -m "Specific topics"
git add ML.md
git commit -m "Machine learning"
git add AI.md
git commit -m "AI"
git add seq.md
git commit -m "sequencing analysis"
for d in DEPICT bowtie2 caviar GATK IGV IMPHLA02 JAM jannovar ldetect ldsc ldpred PASCAL MetaXcan FUSION PLINK2 PyLMM R R-packages PheWAS pindel overall VEP VarScan st.sh
do
   git add $d
   git commit -m "$d"
done
git push
