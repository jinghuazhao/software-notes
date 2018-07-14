# 14-7-2018 JHZ

git add README.md
git commit -m "README"
git add topics.md
git commit -m "Specific topics"
git add ML.md
git commit -m "Machine learning"
git add AI.md
git commit -m "AI"
git add seq.md
git commit -m "Sequencing analysis"
git add misc.md
git commit -m "Miscellaneour programs"
for d in DEPICT caviar IMPHLA02 JAM ldetect ldsc ldpred PASCAL MetaXcan FUSION PLINK2 PyLMM R R-packages PheWAS overall VEP st.sh
do
   git add $d
   git commit -m "$d"
done
git push
