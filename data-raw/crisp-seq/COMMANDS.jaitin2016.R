#!/usr/bin/env Rscript

library(data.table)
library(biomaRt)
library(dplyr)
library(Matrix)
library(R.utils)
library(Rtsne)

redo_from_download <- FALSE
if(redo_from_download) {

    project <- "jaitin2016"
    base_dir <- "~/projects/SCperturb/data-raw/crisp-seq"
    input_dir <- file.path(base_dir, "input")

    # unzip following into input directory
    #https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE90486&format=file
    #UGI_data.zip, personal communication from Assaf Weiner

    output_dir <- "~/projects/SCperturb/data-raw/crisp-seq"

    generate_plots <- FALSE


    # get the ugi to wellid mappings from UGI_data.zip
    ugi_dir <- file.path(input_dir, "UGI")
    ugi_files <- list.files(ugi_dir, pattern=".txt")
    ugi_list <- list()
    for (ugifile in ugi_files){
        ugi_list[[length(ugi_list) + 1]] <- read.table(file.path(ugi_dir, ugifile), header=FALSE, sep=",")
    }
    ugi <- do.call(rbind, ugi_list)
    colnames(ugi) <- c("wellid", "ugi", "umi", "readcount")

    # read the UGI to knock out mapping
    map_file <- file.path(input_dir, "gene_ugi.csv")
    map <- read.table(map_file, header=T, sep="\t")
    colnames(map) <- c("gene", "ugi")
    rownames(map) <- map$UGI

    # keep only rows that are defined in the map
    metadata <- left_join(map, ugi)

    # read the expression matrices
    raw_dir <- file.path(input_dir, "GSE90486_RAW")
    raw_files <- list.files(raw_dir, pattern=".txt.gz")

    col_list <- list()
    expr_list <- list()
    idx <- 1

    for (rawfile in raw_files) {
        d <- read.table(gzfile(file.path(raw_dir, rawfile)), header=T, sep="\t")

        # rows -- reporter genes
        # only need to do this for the first file, because they all have the same rows
        if (idx == 1) {
            ## retireve the gene symbol and chromosomal location
            ensembl <- useMart(host='ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')
            geneInfo <- suppressMessages(getBM(filters = "mgi_symbol", values = rownames(d), attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name"), mart = ensembl))
            # get single symbol for each id
            geneInfo <- as.data.frame(geneInfo %>% group_by(mgi_symbol) %>% summarize(ensembl_gene_id = max(ensembl_gene_id), chromosome_name=max(chromosome_name)))
            rownames(geneInfo) <- geneInfo$mgi_symbol
            rows <- geneInfo[rownames(d),]
            rows$origsymbol <- rownames(d)
            rownames(rows) <- rows$origsymbol
        }


        # columns -- wells
        ugiperwell <- metadata[metadata$wellid %in% colnames(d),]

        # pick the ugi->gene with the most reads
        # count percent of second most abundant reads -- method from rubin2019
        t <- ugiperwell %>% group_by(wellid) %>% summarize(total=sum(readcount)) 
        r <- ugiperwell %>% group_by(wellid) %>% mutate(rank=rank(-readcount))
        j <- left_join(r,t) %>% mutate(perc=readcount/total) 

        # separate files into single perturbations and double perturbations
        # TODO: this method is pretty blunt -- just take the highest number of reads as representative
        single <- as.character(j %>% group_by(wellid) %>% filter(rank==1) %>% pull(wellid))
        cols <- as.data.frame(r %>% filter(wellid %in% single) %>% filter(rank==1) %>% dplyr::select(wellid, readcount, gene) %>% arrange(wellid))
        cols$wellid <- as.character(cols$wellid)
        col_list[[idx]] <- cols
        expr_list[[idx]] <- d[, cols$wellid]
        idx <- idx + 1


        if (generate_plots) {
            ggplot(ugiperwell[], aes(x=ugi)) + geom_histogram(stat="count")
            ggsave(file.path(base_dir, "plots", paste(rawfile, "gg", "pdf", sep=".")))
        
            ugiperwellcounts <- unlist(lapply(colnames(d), function(well) {nrow(ugi[ugi$wellid == well,])}))
            table(ugiperwellcounts)
            pdf(file.path(base_dir, "plots", paste(rawfile, "hist", "pdf", sep=".")))
            hist(ugiperwellcounts, breaks=100)
            dev.off()

            m <- merge(x=ugiperwell, y=map, by="ugi") # TODO: allow errors in barcodes?  , all.x=TRUE)
            m$ugi <- factor(m$ugi, levels=unique(m$ugi))
            ggplot(m, aes(x=readcount, y=gene)) + geom_point(alpha=0.3 ) 
            well <- colnames(d)[[3]]
            m[m$wellid == well,]
            
            u <- m[m$wellid %in% colnames(d),]
            uu <- u %>% group_by(wellid) %>% summarize(count=n())
        }
    }



    c_bind <- do.call(rbind, col_list)
    levels(c_bind$gene) <- sub("control", "CTRL", levels(c_bind$gene))
    d_bind <- do.call(cbind, expr_list)
    mat <- d_bind[,c_bind$wellid]


    write.table(rows, file=file.path(output_dir, paste("rowmetadata", project, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
    write.table(c_bind, file=file.path(output_dir, paste("colmetadata", project, "csv", sep=".")), quote=F, row.names=TRUE, col.names=TRUE, sep="\t")
    sparse <- Matrix(as.matrix(mat), sparse=TRUE)
    sparse_file <- file.path(output_dir, paste("counts", project, "mtx", sep="."))
    writeMM(sparse, file=sparse_file, quote=F, row.names=F, col.names=F)
    if (file.exists(paste(sparse_file, "gz", sep="."))) {
        file.remove(paste(sparse_file, "gz", sep="."))
    }
    gzip(sparse_file)
} else {
    colmetadata.jaitin2016 <- read.table("colmetadata.jaitin2016.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(colmetadata.jaitin2016)

    rowmetadata.jaitin2016 <- read.table("rowmetadata.jaitin2016.csv", sep="\t", quote="", header=TRUE)
    usethis::use_data(rowmetadata.jaitin2016)

    counts.jaitin2016 <- readMM("counts.jaitin2016.mtx.gz")
    usethis::use_data(counts.jaitin2016)
}

