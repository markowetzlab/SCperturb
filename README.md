
# SCperturb
Data and analysis package for single cell perturbation data

# Install and run

Install git-lfs to handle the large files.
Download git-lfs from https://git-lfs.github.com/
or install with brew:

```
$ brew install git-lfs
```

Clone this repository into a local directory, eg ~/projects

Install the SCperturb package locally using devtools:

```
> devtools::install("~/projects/SCperturb")
> library(SCperturb)
```

List available data.  3 data objects are included for each project: counts.project, rowmetadata.project and colmetadata.project.

```
> SCperturb::listprojects()
> data(rowmetadata.datlinger2017_stim)
```

Run mnem on perturb-seq and generate Figure 7 from Holding2019 "Network reconstruction from single cell transcriptomic profiles of CRISPR gene perturbations".  This takes about an hour to run. 

```
> set.seed(42)
> SCperturb::mnem_dixit2016_K562_lowmoi()
```


# Data

For each of the included datasets, this package defines a counts matrix, rowmetadata and colmetdata tables.  Rows describe the reporter genes.  Columms describe the cells, and metadata includes the perturbed genes.  

Data is included for the following papers:

* Dixit2016 (Perturb-seq) 10 perturbations of immune-related transcription factors in K562 cells (GSE90063).

* Datlinger2017 (Crop-seq) 33 perturbations of regulators of TCR signalling and transcription factors in the T-cell receptor pathway in Jurkat cells (GSE92872).

* Jaitin2016 (Crisp-seq) 22 perturbations of immune related genes in mouse (GSE90486)

* Jackson2019 11 perturbations of growth-related genes in 11 different growth conditions in S. cerevisiae (GSE125162) Note that there is only one rowmetadata table for this dataset that applies to all counts matrices.  
