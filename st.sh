# 6-6-2019 JHZ

git add README.md
git commit -m "README"
git add AA.md
git commit -m "Association analysis"
git add CRISPR.md
git commit -m "CRISPR-related topics"
git add SL.md
git commit -m "Statistical learning"
git add AI.md
git commit -m "Artificial intelligence"
git add NGS.md
git commit -m "Next-generation sequencing analysis"
git add prottrans.md
git commit -m "Proteome and transcritome data"
git add misc.md
git commit -m "Miscellaneour programs"
for d in files st.sh
do
   git add $d
   git commit -m "$d"
done
git push
