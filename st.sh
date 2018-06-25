# 25-6-2018 JHZ

git add README.md
git commit -m "README"
git add topics.md
git commit -m "Specific topics"
git add NOTES.md
git commit -m "NOTES"
for d in DEPICT caviar GATK IGV IMPHLA02 JAM jannovar ldsc PASCAL MetaXcan PLINK2 PyLMM R R-packages PheWAS overall VEP VarScan st.sh
do
   git add $d
   git commit -m "$d"
done
git push
