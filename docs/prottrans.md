# Proteome and transcriptome

## CCprofiler/COPF

<https://github.com/CCprofiler/CCprofiler>, <https://sec-explorer.shinyapps.io/hela_cellcycle/>

Heusel M, Bludau I, Rosenberger G, Hafen R, Frank M, Banaei-Esfahani A, Collins B, Gstaiger M, Aebersold R. Complex-centric proteome profiling by SEC-SWATH-MS. *Mol Syst Biol* 2019, 15:e8438 DOI: <https://doi.org/10.15252/msb.20188438>

Heusel et al. A global screen for assembly state changes of the mitotic proteome by SEC-SWATH-MS. *Cell Systems* 0220, 10:1–23, <https://doi.org/10.1016/j.cels.2020.01.001>

## DOGMA

DOGMA is a program for fast and easy quality assessment of transcriptome and proteome data based on conserved protein domains. See

<https://domainworld.uni-muenster.de/programs/dogma/>

<https://ebbgit.uni-muenster.de/domainWorld/DOGMA>

## Immuno-SABER

<https://github.com/HMS-IDAC>

Saka SK, et al. (2019). Immuno-SABER enables highly multiplexed and amplified protein imaging in tissues. *Nat Biotechnol*, <https://www.biorxiv.org/content/10.1101/401810v1>.

## iProFun

<https://github.com/songxiaoyu/iProFun>

Song M, Ji J, Gleason KJ, Yang F, Martignetti JA, Chen LS, Wang P,
Insights into impact of DNA copy number alteration and methylation on the proteogenomic landscape of human ovarian cancer via a multi-omics integrative analysis
*Mol Cell Proteomics*

## MRMassaydb

<http://mrmassaydb.proteincentre.com/fdaassay/>

## Percolator

Semi-supervised learning for peptide identification from shotgun proteomics dataset, and Single Cell Proteomics by Mass Spectrometry (SCOPE-MS).

<http://crux.ms/>, <http://percolator.ms/> and <https://github.com/percolator/percolator>.

McIlwain S, Tamura K, Kertesz-Farkas A, Grant CE, Diament B, Frewen B, Howbert JJ, Hoopmann MR, Käll L, Eng JK, MacCoss MJ, Noble WS (2014). Crux: rapid open source protein tandem mass spectrometry analysis. *J Proteome Res* 13(10):4488-4491.

Fondrie WE, Noble WS (2020). Machine learning strategy that leverages large data sets to boost statistical power in small-scale experiments. *J Proteome Res*

## PeCorA

<https://github.com/jessegmeyerlab/PeCorA>

## PIVar

<https://github.com/WeiWenqing/PIVar>

Teng H, et al. Prevalence and architecture of posttranscriptionally impaired synonymous mutations in 8,320 genomes across 22 cancer types [published online ahead of print, 2020 Jan 17]. Nucleic Acids Res. 2020;gkaa019. doi:10.1093/nar/gkaa019

##

[PROPERseqTools](https://github.com/Zhong-Lab-UCSD/PROPERseqTools/tree/v1.0.0) and [database](https://genemo.ucsd.edu/proper/)

Johnson, K.L. et al. Revealing protein-protein interactions at the transcriptome scale by sequencing. Molecular Cell.

## ProteomeXchange

<http://proteomecentral.proteomexchange.org/cgi/GetDataset>, also PRotein IDEntification databases (PRIDE): <https://www.ebi.ac.uk/>

## ProteoformAnalysis

<https://github.com/ibludau/ProteoformAnanlysis>

<http://proteoformviewer.ethz.ch/>, <http://proteoformviewer.ethz.ch/ProteoformExplorer_final_classes>

```bash
# In silico benchmark.
install_github("CCprofiler/CCprofiler", ref =  "proteoformLocationMapping")
R --no-save -q < InterlabBenchmark_final_paper.R
# COPF analysis of the cell cycle SEC-SWATH-MS dataset
## E1709051521_feature_alignment.tsv
wget -qO- http://ftp.pride.ebi.ac.uk/pride/data/archive/2019/05/PXD010288/SWATH_QueryResultGlobalAlignment_20170906213658025-1332538.zip | \
unzip SWATH_QueryResultGlobalAlignment_20170906213658025-1332538.zip
wget http://ftp.pride.ebi.ac.uk/pride/data/archive/2019/05/PXD010288/HeLaCCL2_SEC_annotation_full.xlsx
...
```

Bludau I. et al. Systematic detection of functional proteoform groups from bottom-up proteomic datasets. *Nat Comm* 2021, 12:3810

## QTLtools

<https://qtltools.github.io/qtltools/>

## RNAsnp

<https://rth.dk/resources/rnasnp/software>

## ViennaRNA

<https://www.tbi.univie.ac.at/RNA/>
