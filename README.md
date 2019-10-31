
# SCperturb
Data and analysis package for single cell perturbation data to complement Holding2019 "Network reconstruction from single cell transcriptomic profiles of CRISPR gene perturbations".  This package includes data in a standardized format from 4 single cell perturbation experiments, and example code to construct a network from Datlinger2017 using mnem (Pirkl2018).  


# Install and run

An R markdown file, Holding2019.Rmd, is also available that covers the following steps. 

Install the SCperturb package from github using devtools.  When the package is installed, it will download 285Mb of additional data from https://content.cruk.cam.ac.uk/fmlab/holding2019

Note: this package depends on the pcalg package, which requires that a fortran compiler be installed, and installation will give confusing error messages if this requirement isn't met.  

```
> install.packages(devtools)
> devtools::install_github("markowetzlab/SCperturb")
> library(SCperturb)
```

List available data.  3 data objects are included for each project: counts.project, rowmetadata.project and colmetadata.project.

```
> SCperturb::listprojects()
> data(rowmetadata.datlinger2017_stim)
```

Run mnem on perturb-seq data and generate Figure 7 from Holding2019.  This takes about an hour to run. 

```
> set.seed(42)
> SCperturb::mnem_dixit2016_K562_lowmoi()
```
Generate Figure 6 from Holding2019 that shows the knockdown efficiency.

```
> knockdown_efficiency()
```

# Data

For each of the included datasets, this package defines a counts matrix, rowmetadata and colmetdata tables.  Rows describe the reporter genes.  Columms describe the cells, and metadata includes the perturbed genes.  

Data is included for the following papers:

* Dixit2016 (Perturb-seq) 10 perturbations of immune-related transcription factors in K562 cells (GSE90063).

* Datlinger2017 (Crop-seq) 33 perturbations of regulators of TCR signalling and transcription factors in the T-cell receptor pathway in Jurkat cells (GSE92872).

* Jaitin2016 (Crisp-seq) 22 perturbations of immune related genes in mouse (GSE90486)

* Jackson2019 11 perturbations of growth-related genes in 11 different growth conditions in S. cerevisiae (GSE125162) Note that there is only one rowmetadata table for this dataset that applies to all counts matrices.  



# Regenerate R data objects

If you want to reproduce the R data objects that are distributed with this package from scratch, then you can find these scripts in /data-raw in the package source.  Large source files can be downloaded from https://content.cruk.cam.ac.uk/fmlab/holding2019/data-raw/.
