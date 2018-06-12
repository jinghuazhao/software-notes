# 12-6-2018 JHZ

git add README.md
git commit -m "README"
for d in DEPICT caviar JAM jannovar ldsc PASCAL MetaXcan PLINK2 PyLMM R R-packages PheWAS overall VEP VarScan st.sh
do
   git add $d
   git commit -m "$d"
done
git push
