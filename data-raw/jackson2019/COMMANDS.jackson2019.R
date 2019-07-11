
library(dplyr)
library(biomaRt)
library(Matrix)
library(R.utils)
library(data.table)

project <- "jackson2019"
conditions <- c("YPD", "DIAUXY", "RAPA", "YPETOH", "MMD", "MMETOH", "NLIMGLN", "NLIMPRO", "NLIMNH4", "NLIMUREA", "CSTARVE")
redo_from_download <- FALSE

if(redo_from_download) {
    base_dir <- "~/projects/SCperturb/data-raw/jackson2019"
    input_dir <- file.path(base_dir, "input")

    # untar the following files into the input dir
    # https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE125162&format=file

    ensembl <- useMart(host='ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='scerevisiae_gene_ensembl')

    setwd(input_dir)
    needrows <- TRUE
    for (cond in conditions) {
        metadata_file <- list.files(pattern=paste(cond, "-fastqTomat0-BCLINK.tsv.gz", sep=""))[[1]]
        metadata <- read.table(gzfile(file.path(input_dir, metadata_file)))
        counts_file <- list.files(pattern=paste(cond, "-fastqTomat0-Counts.tsv.gz", sep=""))[[1]]
        counts <- read.table(gzfile(file.path(input_dir, counts_file)))
        ct <- counts[,-(6830:6832)]
        sparse <- Matrix(t(as.matrix(ct)))


        rows <- data.frame()
        if (needrows) {
            # rows -- reporter genes
            rows <- data.frame(ensembl_gene_id=rownames(sparse))
            geneInfo <- getBM(filters = "ensembl_gene_id", values = rownames(sparse), attributes = c("ensembl_gene_id", "chromosome_name"), mart = ensembl)
            # get single symbol for each id
            geneInfo <- as.data.frame(geneInfo %>% group_by(ensembl_gene_id) %>% summarize(chromosome_name=max(chromosome_name)))
            rows <- left_join(rows, geneInfo)
            rownames(rows) <- rows$ensembl_gene_id
            needrows <- FALSE
            write.table(rows, file=file.path(base_dir, paste("rowmetadata", project, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
        }
        

        # columns - cells
        cols <- metadata[metadata$Cell_Barcode %in% colnames(sparse),]
        cols <- metadata[metadata$Cell_Barcode %in% rownames(ct),]
        # check these have unique genotypes
        n <- cols %>% group_by(Cell_Barcode) %>% summarize(n=n_distinct(Genotype))
        stopifnot(all.equal(which(n$n != 1), integer(0)))
        # all the barcodes that have entries in the counts matrix describe single perturbations, but may be listed multiple times, so take sum
        cols <- as.data.frame(cols %>% group_by(Cell_Barcode) %>% summarize(Genotype=max(as.character(Genotype)), UMI_Count=sum(UMI_Count)) %>% rename(grna = Genotype))
        cols$gene <- as.factor(unlist(lapply(as.character(cols$grna), function(x) {strsplit(x, "_")[[1]][[1]]})))
        rownames(cols) <- cols$Cell_Barcode


        thisproject <- paste(project, cond, sep="_")

        write.table(cols, file=file.path(base_dir, paste("colmetadata", thisproject, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
        sparse_file <- file.path(base_dir, paste("counts", thisproject, "mtx", sep="."))
        writeMM(sparse, file=sparse_file, quote=F, row.names=F, col.names=F)
        if (file.exists(paste(sparse_file, "gz", sep="."))) {
            file.remove(paste(sparse_file, "gz", sep="."))
        }
        gzip(sparse_file)
    }
} else {
    rowmetadata.jackson2019 <- read.table("rowmetadata.jackson2019.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(rowmetadata.jackson2019)

    colmetadata.jackson2019_YPD <- read.table("colmetadata.jackson2019_YPD.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_YPD)
    colmetadata.jackson2019_DIAUXY <- read.table("colmetadata.jackson2019_DIAUXY.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_DIAUXY)
    colmetadata.jackson2019_RAPA <- read.table("colmetadata.jackson2019_RAPA.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_RAPA)
    colmetadata.jackson2019_YPETOH <- read.table("colmetadata.jackson2019_YPETOH.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_YPETOH)
    colmetadata.jackson2019_MMD <- read.table("colmetadata.jackson2019_MMD.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_MMD)
    colmetadata.jackson2019_MMETOH <- read.table("colmetadata.jackson2019_MMETOH.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_MMETOH)
    colmetadata.jackson2019_NLIMGLN <- read.table("colmetadata.jackson2019_NLIMGLN.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_NLIMGLN)
    colmetadata.jackson2019_NLIMPRO <- read.table("colmetadata.jackson2019_NLIMPRO.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_NLIMPRO)
    colmetadata.jackson2019_NLIMNH4 <- read.table("colmetadata.jackson2019_NLIMNH4.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_NLIMNH4)
    colmetadata.jackson2019_NLIMUREA <- read.table("colmetadata.jackson2019_NLIMUREA.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_NLIMUREA)
    colmetadata.jackson2019_CSTARVE <- read.table("colmetadata.jackson2019_CSTARVE.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jackson2019_CSTARVE)

    counts.jackson2019_YPD <- readMM("counts.jackson2019_YPD.mtx.gz")
    usethis::use_data(counts.jackson2019_YPD)
    counts.jackson2019_DIAUXY <- readMM("counts.jackson2019_DIAUXY.mtx.gz")
    usethis::use_data(counts.jackson2019_DIAUXY)
    counts.jackson2019_RAPA <- readMM("counts.jackson2019_RAPA.mtx.gz")
    usethis::use_data(counts.jackson2019_RAPA)
    counts.jackson2019_YPETOH <- readMM("counts.jackson2019_YPETOH.mtx.gz")
    usethis::use_data(counts.jackson2019_YPETOH)
    counts.jackson2019_MMD <- readMM("counts.jackson2019_MMD.mtx.gz", sep="\t", quote="", header=TRUE)
    usethis::use_data(counts.jackson2019_MMD)
    counts.jackson2019_MMETOH <- readMM("counts.jackson2019_MMETOH.mtx.gz", sep="\t", quote="", header=TRUE)
    usethis::use_data(counts.jackson2019_MMETOH)
    counts.jackson2019_NLIMGLN <- readMM("counts.jackson2019_NLIMGLN.mtx.gz", sep="\t", quote="", header=TRUE)
    usethis::use_data(counts.jackson2019_NLIMGLN)
    counts.jackson2019_NLIMPRO <- readMM("counts.jackson2019_NLIMPRO.mtx.gz", sep="\t", quote="", header=TRUE)
    usethis::use_data(counts.jackson2019_NLIMPRO)
    counts.jackson2019_NLIMNH4 <- readMM("counts.jackson2019_NLIMNH4.mtx.gz", sep="\t", quote="", header=TRUE)
    usethis::use_data(counts.jackson2019_NLIMNH4)
    counts.jackson2019_NLIMUREA <- readMM("counts.jackson2019_NLIMUREA.mtx.gz", sep="\t", quote="", header=TRUE)
    usethis::use_data(counts.jackson2019_NLIMUREA)
    counts.jackson2019_CSTARVE <- readMM("counts.jackson2019_CSTARVE.mtx.gz", sep="\t", quote="", header=TRUE)
    usethis::use_data(counts.jackson2019_CSTARVE)
}

