#!/usr/bin/bash

function old()
{
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
  git add pharmacogenomics.md
  git commit -m "Pharmacogenomics"
  git add prottrans.md
  git commit -m "Proteome and transcritome data"
  git add misc.md
  git commit -m "Miscellaneour programs"
  git add scRNASeq.md
  git commit -m "Single Cell RNA-Seq"
  git add bin
  git commit -m "PDF editor 1.39, quanto 1.2.4, (un)rar 5.50"
  git add files
  git commit -m "files"
  git push
}

function setup()
{
  module load python/3.8
  source ~/rds/public_databases/software/py38/bin/activate
}

if [ "$(uname -n | sed 's/-[0-9]*$//')" == "login-q" ]; then
   module load ceuadmin/libssh/0.10.6-icelake
   module load ceuadmin/openssh/9.7p1-icelake
   module load ceuadmin/krb5
fi

setup
mkdocs build
mkdocs gh-deploy
git add .gitignore
git add docs
git add mkdocs.yml
git commit -m "backup"
git push
