# 16-7-2018 JHZ

git add README.md
git commit -m "README"
git add association.md
git commit -m "Specific topics"
git add ML.md
git commit -m "Machine learning"
git add AI.md
git commit -m "AI"
git add seq.md
git commit -m "Sequencing analysis"
git add misc.md
git commit -m "Miscellaneour programs"
for d in MetaXcan R-packages envirs st.sh
do
   git add $d
   git commit -m "$d"
done
git push
