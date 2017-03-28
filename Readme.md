# R code used for data analysis of Jimenez-Martinez *et al.*, Science (2017)

The R scripts for read mapping and counting, computational filtering (quality control), normalisation, differential expression analysis and for figure generation can be found in this repository. 

* The Makefile for read alignment and counting can be found in [Mapping](https://github.com/MarioniLab/ImmuneAging2017/tree/master/Mapping/).

* Quality control and filtering of low quality cells/genes is described in [Quality_control](https://github.com/MarioniLab/ImmuneAging2017/tree/master/Quality_control/). This folder contains the scripts to reproduce Fig. S1 and Fig. S2.

* [Normalisation](../master/Normalisation/) contains scripts to normalise scRNAseq data using the [BASiCS](https://github.com/catavallejos/BASiCS) package.

* The [Differential_testing](https://github.com/MarioniLab/ImmuneAging2017/tree/master/Differential_testing/) repository holds scripts that perform differential expression and differential variability testing of scRNAseq data.

* All scripts to reproduce the main and supplementary figures can be found in [Analysis](https://github.com/MarioniLab/ImmuneAging2017/tree/master/Analysis/)

* Raw and normalised transcript counts as well as some metadata information are available in [Data](https://github.com/MarioniLab/ImmuneAging2017/tree/master/Data/)
