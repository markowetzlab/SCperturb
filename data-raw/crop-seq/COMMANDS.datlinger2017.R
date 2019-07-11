
library(Matrix)
library(biomaRt)
library(dplyr)
library(R.utils)

redo_from_download <- FALSE
if(redo_from_download) {

    project <- "datlinger2017"
    base_dir <- "~/projects/SCperturb/data-raw/crop-seq"
    input_dir <- file.path(base_dir, "input")

    # put following file into input dir
    # https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE92872&format=file&file=GSE92872%5FCROP%2Dseq%5FJurkat%5FTCR%2Edigital%5Fexpression%2Ecsv%2Egz

    filename <- "GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv.gz"
    # read same file to extract metadata then get the data
    metadata <- read.table(gzfile(file.path(input_dir, filename)), header=F, sep=",", nrows=5)
    data <- read.table(gzfile(file.path(input_dir, filename)), header=F, sep=",", skip=6)
    rownames(data) <- data$V1
    data <- data[,-1]


    # columns == cells

    cells <- as.data.frame(t(metadata)[-1,])
    colnames(cells) <- t(metadata)[1,]
    cells$condition <- factor(cells$condition)
    # get rid of cells that didn't read uniquely
    exclude <- which(cells$cell == "GTCGTTTTAACN")
    cells <- cells[-exclude,]
    data <- data[,-exclude]
    rownames(cells) <- cells$cell

    cond <- which(cells$condition == "stimulated")
    cells_stim <- cells[cond,]
    cells_unstim <- cells[-cond,]

    data_stim <- data[,cond]
    data_unstim <- data[,-cond]

    # rows == reporter genes

    # retireve the gene symbol and chromosomal location
    ensembl <- useMart(host='ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
    geneInfo <- suppressMessages(getBM(filters = "hgnc_symbol", values = rownames(data), attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"), mart = ensembl))

    # the data is unfortunately indexed by gene symbol
    rows <- as.data.frame(geneInfo %>% group_by(hgnc_symbol) %>% summarize(ensembl_gene_id=max(ensembl_gene_id), chromosome_name=max(chromosome_name)))
    rownames(rows) <- rows$hgnc_symbol
    # get single symbol for each id
    rows <- rows[rownames(data),]
    rows$origsymbol <- rownames(data)
    rownames(rows) <- rownames(data)


    proj <- paste(project, "stim", sep="_")
    write.table(rows, file=file.path(base_dir, paste("rowmetadata", proj, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
    write.table(cells_stim, file=file.path(base_dir, paste("colmetadata", proj, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
    sparse <- Matrix(as.matrix(data_stim), sparse=TRUE)
    sparse_file <- file.path(base_dir, paste("counts", proj, "mtx", sep="."))
    writeMM(sparse, file=sparse_file, quote=F, row.names=FALSE, col.names=FALSE)
    if (file.exists(paste(sparse_file, "gz", sep="."))) {
        file.remove(paste(sparse_file, "gz", sep="."))
    }
    gzip(sparse_file)

    proj <- paste(project, "unstim", sep="_")
    write.table(rows, file=file.path(base_dir, paste("rowmetadata", proj, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
    write.table(cells_unstim, file=file.path(base_dir, paste("colmetadata", proj, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
    sparse <- Matrix(as.matrix(data_unstim), sparse=TRUE)
    sparse_file <- file.path(base_dir, paste("counts", proj, "mtx", sep="."))
    writeMM(sparse, file=sparse_file, quote=F, row.names=FALSE, col.names=FALSE)
    if (file.exists(paste(sparse_file, "gz", sep="."))) {
        file.remove(paste(sparse_file, "gz", sep="."))
    }
    gzip(sparse_file)
    #sce <- SingleCellExperiment(assays = list(counts=sparse), colData=cells, rowData=rows)
}
else {
    colmetadata.datlinger2017_stim <- read.table("colmetadata.datlinger2017_stim.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.datlinger2017_stim)

    rowmetadata.datlinger2017_stim <- read.table("rowmetadata.datlinger2017_stim.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(rowmetadata.datlinger2017_stim)

    counts.datlinger2017_stim <- readMM("counts.datlinger2017_stim.mtx.gz")
    usethis::use_data(counts.datlinger2017_stim)

    colmetadata.datlinger2017_unstim <- read.table("colmetadata.datlinger2017_unstim.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.datlinger2017_unstim)

    rowmetadata.datlinger2017_unstim <- read.table("rowmetadata.datlinger2017_unstim.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(rowmetadata.datlinger2017_unstim)

    counts.datlinger2017_unstim <- readMM("counts.datlinger2017_unstim.mtx.gz")
    usethis::use_data(counts.datlinger2017_unstim)
}

