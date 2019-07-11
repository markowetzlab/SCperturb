
library(data.table)
library(reshape2)
library(dplyr)
library(Matrix)
library(biomaRt)
library(R.utils)


# cols -- functions

readmap <- function(input_file) {
    # use gbc to cbc mapping to fix resolve the multiple perturbations
    map <- read.table(input_file, header=FALSE, quote="\"", sep=",", stringsAsFactors=FALSE)
    longmap <- list()
    longmapidx <- 1
    for (rowidx in 1:nrow(map)) {
        y <- strsplit(map[rowidx,2], ", ")[[1]]
        longmap <- append(longmap, lapply(y, function(x) {c(map[rowidx,1], x)}))
    } 
    map <- as.data.frame(do.call(rbind, longmap))
    colnames(map) <- c("guide", "cbc")
    #rownames(map) <- map$cbc
    return(map)
}

binary_exp_mat <- function() {
    #construct Binary experiment matrix. Each row corresponds to an cell/experiment, each column to a pathway gene.  The entry in row e and column g is equal 1 if experiment e targets gene g.
    map <- readmap(file.path(input_dir, "GSM2396860_k562_tfs_highmoi_cbc_gbc_dict_lenient.csv.gz"))
    wide <- dcast(map, cbc ~ guide)
    cbcs<- wide$cbc
    wide$cbc <- NULL
    wide[is.na(wide)] <- 0
    wide[wide != 0] <- 1
    wide <- sapply(wide, as.numeric)
    rownames(wide) <- cbcs
    #wide[1:6,1:6]
    wide <- as.data.frame(wide)

    # aggregate columns from guides to genes
    guides <- colnames(wide)
    guides <- gsub("_[0-9]+", "", guides)
    guides <- gsub("p_sg", "", guides)

    guides <- gsub("p_IN", "", guides)
    guides <- gsub("NIC", "", guides)
    guides <- gsub("TERGE", "CTRL", guides)

    colnames(wide) <- guides
    aggr <- by(t(wide), INDICES=row.names(t(wide)), FUN=colSums)
    wide <- as.data.frame(do.call(cbind,aggr))

    sgenes <- colnames(wide)[5:14]
    w <- wide[5:14]
    w[w>1] <- 1


    # do we have at least one experiment (cell) per single gene?
    found1 <- 0
    for (c1 in sgenes) {
        if (length(which(w[,c1] == 1 & unname(rowSums(w) == 1))) >= 1) {
            found1 <- found1 + 1
        } else {
            print(c1, "single perturbation missing", sep=" ")
        }
    }
    print(paste("singles:", found1, "of 10 present", sep=" "))
    length(which(unname(rowSums(wide) == 1)))

    # do we have at least one cell per double perturbation?
    found2 <- 0
    c <- combn(sgenes, 2)
    for (idx in 1:ncol(c)) {
        c1 <- c[1, idx]
        c2 <- c[2, idx]
        if (length(which(w[[c1]] == 1 & w[[c2]] == 1 & unname(rowSums(w) == 2))) >= 1) {
            found2 <- found2 + 1
        } else {
            print(paste(c1, c2, "double perturbation missing", sep=" ")) 
        }
    }
    print(paste("doubles:", found2, "of", ncol(c), "present"))
    length(which(unname(rowSums(wide) == 2)))


    found3 <- 0
    c <- combn(sgenes, 3)
    for (idx in 1:ncol(c)) {
        c1 <- c[1, idx]
        c2 <- c[2, idx]
        c3 <- c[3, idx]
        if (length(which(w[[c1]] == 1 & w[[c2]] == 1 & w[[c3]] == 1 & unname(rowSums(w) == 3))) >= 1) {
            found3 <- found3 + 1
        } else {
            print(paste(c1, c2, c3, "triple perturbation missing", sep=" ")) 
        }
    }
    print(paste("triples:", found3, "of", ncol(c), "present", sep=" "))
    length(which(unname(rowSums(wide) == 3)))

    found4 <- 0
    c <- combn(sgenes, 4)
    for (idx in 1:ncol(c)) {
        c1 <- c[1, idx]
        c2 <- c[2, idx]
        c3 <- c[3, idx]
        c4 <- c[4, idx]
        if (length(which(w[[c1]] == 1 & w[[c2]] == 1 & w[[c3]] == 1 & w[[c4]] == 1 & unname(rowSums(w) == 4))) >= 1) {
            found4 <- found4 + 1
        }
    }
    print(paste("quads: ", found4, "of", ncol(c), "present", sep=" "))
    length(which(unname(rowSums(wide) == 4)))

    found5 <- 0
    c <- combn(sgenes, 5)
    for (idx in 1:ncol(c)) {
        c1 <- c[1, idx]
        c2 <- c[2, idx]
        c3 <- c[3, idx]
        c4 <- c[4, idx]
        c5 <- c[5, idx]
        if (length(which(w[[c1]] == 1 & w[[c2]] == 1 & w[[c3]] == 1 & w[[c4]] == 1 & w[[c5]] == 1 & unname(rowSums(w) == 5))) >= 1) {
            found5 <- found5 + 1
        }
    }
    print(paste("pents: ", found5, "of", ncol(c), "present", sep=" "))
    length(which(unname(rowSums(wide) == 5)))


    return(wide)
}


redo_from_download <- FALSE
if(redo_from_download) {

    project <- "dixit2016_K562_highmoi"
    base_dir <- "~/projects/sc_perturb/data-raw/perturb-seq"
    input_dir <- file.path(base_dir, "input")

    # untar the following files into the input directory
    # https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE90063&format=file


    sparse <- readMM(file.path(input_dir, "GSM2396860_k562_tfs_highmoi.mtx.txt.gz"))
    genes <- read.table(file.path(input_dir, "GSM2396860_k562_tfs_highmoi_genenames.csv.gz"), sep=",", header=TRUE)
    cells <- read.table(file.path(input_dir, "GSM2396860_k562_tfs_highmoi_cellnames.csv.gz"), sep=",", header=TRUE)
    stopifnot(all.equal(nrow(sparse), nrow(genes)))
    stopifnot(all.equal(ncol(sparse), nrow(cells)))

    stopifnot(all.equal(cells[,1] + 1, 1:nrow(cells)))
    stopifnot(all.equal(genes[,1] + 1, 1:nrow(genes)))
    rownames(cells) <- cells[,2]
    rownames(genes) <- genes[,2]


    # rows -- reporter genes

    genes$ensg <- unlist(lapply(rownames(genes), function(x) {strsplit(x, "_")[[1]][[1]]}))
    rownames(genes) <- genes$ensg
    ## retireve the gene symbol and chromosomal location
    ensembl <- useMart(host='useast.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
    geneInfo <- suppressMessages(getBM(filters = "ensembl_gene_id", values = genes$ensg, attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"), mart = ensembl))
    # get single symbol for each id
    geneInfo <- as.data.frame(geneInfo %>% group_by(ensembl_gene_id) %>% summarize(hgnc_symbol = max(hgnc_symbol), chromosome_name=max(chromosome_name)))
    rows <- left_join(genes, geneInfo, by=c("ensg"="ensembl_gene_id"))
    rownames(rows) <- rows$ensg


    # cols -- cells

    # get the genes that are perturbed in each cell
    cols <- binary_exp_mat()
    # put them in the right order
    cols <- cols[rownames(cells),]
    
    #why does one row end up all na?  we don't want it anyway
    cols <- cols[-321,]
    sparse <- sparse[,-321]

    colnames(cols) <- paste("gene", colnames(cols), sep="_")
    cols$cell <- rownames(cols)


    write.table(rows, file=file.path(base_dir, paste("rowmetadata", project, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
    write.table(cols, file=file.path(base_dir, paste("colmetadata", project, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
    sparse_file <- file.path(base_dir, paste("counts", project, "mtx", sep="."))
    writeMM(sparse, file=sparse_file, quote=F, row.names=F, col.names=F)
    if (file.exists(paste(sparse_file, "gz", sep="."))) {
        file.remove(paste(sparse_file, "gz", sep="."))
    }
    gzip(sparse_file)
}
else {
    colmetadata.dixit2016_K562_highmoi <- read.table("colmetadata.dixit2016_K562_highmoi.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.dixit2016_K562_highmoi)

    rowmetadata.dixit2016_K562_highmoi <- read.table("rowmetadata.dixit2016_K562_highmoi.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(rowmetadata.dixit2016_K562_highmoi)

    counts.dixit2016_K562_highmoi <- readMM("counts.dixit2016_K562_highmoi.mtx.gz")
    usethis::use_data(counts.dixit2016_K562_highmoi)
}


