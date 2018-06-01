## Issues with DEPICT

## Setup

* Information on download and usage is here, https://data.broadinstitute.org/mpg/depict/documentation.html, which has links on [DEPICT_v1_rel194.tar.gz](https://data.broadinstitute.org/mpg/depict/depict_download/bundles/DEPICT_v1_rel194.tar.gz).
It contains 1000Genomes and other data unavailable from [depict_140721.tar.bz2](https://data.broadinstitute.org/mpg/depict/depict_140721.tar.bz2).
* It is preferable to use the source package from GitHub while linking the data directory packaged with DEPICT_v1_rel194.tar.gz above. For instance, 
```{bash}
tar xvfz DEPICT_v1_rel194.tar.gz
export CWD=$(pwd)
git clone https://github.com/perslab/depict
cd depict
mv data data.sav
ln -s $CWD/DEPICT/data
ln -s src/python/depict.py $HOME/bin/depict.py
cd example
# editing ldl_teslovich_nature2010.cfg
depict.py ldl_teslovich_nature2010.cfg
```
where the package is unpacked into the DEPICT/ directory which contains data/ subdirectory. We then create a symbolic link to data/ from the GitHub version, which does not contain data (to save space on GitHub) while allowing for cutoff_type to be p-values in network analysis for instance.

Note that the documentation example does not give the full results; in order to do so change the following section,
```
# The reconstituted gene set files used by DEPICT
reconstituted_genesets_file: data/reconstituted_genesets/reconstituted_genesets_example.txt
```
to
```
# The reconstituted gene set files used by DEPICT
reconstituted_genesets_file: data/reconstituted_genesets/reconstituted_genesets_150901.binary
```

Under Windows, `gzip.exe` is also required at the working directory or %path%. We can then execute
```
python depict.py BMI.cfg
```

For tissue plot, one can use pdftopng from XpdfReader (or convert/magick from ImageMagick) to obtain .png files to be incorporated into Excel workbook. For network plot, the python package scikit-learn is required.
```bash
sudo pip install scikit-learn
```

## PLINK

[PLINK-1.9](https://www.cog-genomics.org/plink2/), with --clump option, has to be used -- [PLINK2](https://www.cog-genomics.org/plink/2.0/) drops the --clump option.

## Python 2.7.*

For install, the following change is needed: from .sort() to .sort_values() in network_plot.py and depict_library.py .

Is it possible to replicate the Supplementary Figure 9 of the Scott paper? The number of significant pathways seemed to fall short of the FDR<=0.05 criterion. See
[SUMSTATS](https://github.com/jinghuazhao/SUMSTATS) for how to set up.

## template.cfg

This is from src/python rather than .cfg from example.

## Additional notes

You can examine my [PW-pipeline](https://github.com/jinghuazhao/PW-pipeline) repository on other changes I have made.

## Example

We illustrate with the latest GIANT+Biiobank data on BMI,

```bash
wget https://portals.broadinstitute.org/collaboration/giant/images/0/0f/Meta-analysis_Locke_et_al%2BUKBiobank_2018.txt.gz
gunzip -c Meta-analysis_Locke_et_al+UKBiobank_2018.txt.gz | awk '
{
   FS=OFS="\t"
   if(NR==1) print "SNP","Chr","Pos","P"
   else print $3,$1,$2,$9
}' | gzip -f > BMI.txt.gz

```
where we opt to customise the header rather than the DEPICT configuration file.

Once started, we had complaint that
```
Retrieving background loci
Exiting.. To few background files in data/backgrounds/nloci723_nperm500_kb500_rsq0.1_mhc25000000-35000000_colld0.5-collection-1000genomespilot-depict-150429/. Please remove the folder, rerun DEPICT and contact tunepers@broadinstitute.org if the error prevails.
```
Follow instruction and remove the directory. It is very slow-going, ~20 hours on our Linux note but surprisingly half that time under my Windows 10 whose directory zipped and then unzipped under Linux and run`depict.py`there.

We then generate [BMI.xlsx](BMI.xlsx) as in [PW-pipelne](https://github.com/jinghuazhao/PW-pipeline). While there are 849 genesets with FDR<0.05, tissue enrichment shows compelingly an overwhelming role of the nervous system.
