#!/usr/bin/env Rscript

# this is the first, low moi experiment

library(data.table)
library(biomaRt)
library(dplyr)
library(Matrix)
library(R.utils)

redo_from_download <- FALSE
if(redo_from_download) {

    project <- "dixit2016_K562_lowmoi"
    base_dir <- "~/projects/SCperturb/data-raw/perturb-seq"
    input_dir <- file.path(base_dir, "input")

    # put the following files into the input dir
    # k652_both_filt.txt
    # k652_metadata.txt
    # download from Broad Single Cell Portal, Study: Perturb-seq K562
    # https://portals.broadinstitute.org/single_cell/study/SCP31/perturb-seq-k562#study-download

    expression <- as.data.frame(fread(file=file.path(input_dir, "k562_both_filt.txt"), header=TRUE, stringsAsFactors=FALSE))
    rownames(expression) <- expression$GENE
    expression <- expression[,-1] # leave off the gene column


    # columns == cells

    meta <- read.table(file.path(input_dir, "k562_metadata.txt"), header=FALSE, skip=2)
    colnames(meta) <- c("cell", "pool", "perturbation")
    meta$gene <- gsub("_[0-9]+", "", meta$perturbation)
    meta$gene <- gsub("sg", "", meta$gene)
    meta$gene <- gsub("INTERGENIC", "CTRL", meta$gene)
    rownames(meta) <- meta$cell
    # keep only promoter experiments
    meta_prom <- meta[meta$pool == "promoters",]
    meta_prom$cell <- factor(meta_prom$cell)
    meta_prom$pool <- factor(meta_prom$pool)
    meta_prom$perturbation <- factor(meta_prom$perturbation)
    # keep only single perturbations
    meta_prom <- meta_prom[meta_prom$perturbation != "multi",]


    # rows == reporter genes

    # retireve the gene symbol and chromosomal location
    ensembl <- useMart(host='ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
    geneInfo <- suppressMessages(getBM(filters = "hgnc_symbol", values = rownames(expression), attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"), mart = ensembl))

    # the data is unfortunately indexed by gene symbol
    rows <- as.data.frame(geneInfo %>% group_by(hgnc_symbol) %>% summarize(ensembl_gene_id=max(ensembl_gene_id), chromosome_name=max(chromosome_name)))
    rownames(rows) <- rows$hgnc_symbol
    rows <- rows[rownames(expression),]
    rows$origsymbol <- rownames(expression)
    rownames(rows) <- rownames(expression)

    # select expression data based on the rows and columns we have specified above
    expression_tf <- expression[rownames(rows), rownames(meta_prom)]


    write.table(rows, file=file.path(base_dir, paste("rowmetadata", project, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
    write.table(meta_prom, file=file.path(base_dir, paste("colmetadata", project, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
    sparse <- Matrix(as.matrix(expression_tf), sparse=TRUE)
    rownames(sparse) <- rownames(rows)
    sparse_file <- file.path(base_dir, paste("counts", project, "mtx", sep="."))
    writeMM(sparse, file=sparse_file, quote=F, row.names=F, col.names=F)
    if (file.exists(paste(sparse_file, "gz", sep="."))) {
        file.remove(paste(sparse_file, "gz", sep="."))
    }
    gzip(sparse_file)
}
else {
    colmetadata.dixit2016_K562_lowmoi <- read.table("colmetadata.dixit2016_K562_lowmoi.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.dixit2016_K562_lowmoi)

    rowmetadata.dixit2016_K562_lowmoi <- read.table("rowmetadata.dixit2016_K562_lowmoi.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(rowmetadata.dixit2016_K562_lowmoi)

    counts.dixit2016_K562_lowmoi <- readMM("counts.dixit2016_K562_lowmoi.mtx.gz")
    usethis::use_data(counts.dixit2016_K562_lowmoi)

}

#library(SingleCellExperiment)
#ps <- SingleCellExperiment(assays = list(counts = as.matrix(expression_tf)), colData=meta_prom, rowData=rows)

