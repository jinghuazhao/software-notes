# 16-7-2018 JHZ

git add README.md
git commit -m "README"
git add association.md
git commit -m "Specific topics"
git add SL.md
git commit -m "Statistical learning"
git add AI.md
git commit -m "AI"
git add NGS.md
git commit -m "Next-generation sequencing analysis"
git add misc.md
git commit -m "Miscellaneour programs"
for d in envirs st.sh
do
   git add $d
   git commit -m "$d"
done
git push
