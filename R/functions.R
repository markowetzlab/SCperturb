
#!/usr/bin/env Rscript


#' Generate figure showing perturb-seq K562 lowmoi network constructed with mnem
#'
#' @return nothing
#' @export
mnem_dixit2016_K562_lowmoi <- function() {
    projects <- listprojects()
    project <- projects[[4]]

    # get data as data objects from package
    data(colmetadata.dixit2016_K562_lowmoi)
    data(rowmetadata.dixit2016_K562_lowmoi)
    data(counts.dixit2016_K562_lowmoi)
    colmetadata <- colmetadata.dixit2016_K562_lowmoi
    rowmetadata <- rowmetadata.dixit2016_K562_lowmoi
    countsdata <- counts.dixit2016_K562_lowmoi

    # do quality control
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=as.matrix(countsdata)), colData=colmetadata, rowData=rowmetadata)
    sce <- qualitycontrol(sce)
    countsdata <- DESeq2::counts(sce)
    colmetadata <- SummarizedExperiment::colData(sce)
    rowmetadata <- SummarizedExperiment::rowData(sce)

    # calculate differential expression
    p_threshold <- 0.00000005
    diffexp <- diffexp_seurat(countsdata, rowmetadata, colmetadata)
    filtered <- diffexp[diffexp$p_val_adj < p_threshold,]
    diffexp_genes <- unique(filtered$egene)

    # format data for input to mnem
    selected_rows <- which(rowmetadata$origsymbol %in% diffexp_genes)
    selected_counts <- countsdata[selected_rows,]
    rownames(selected_counts) <- rowmetadata[selected_rows,]$origsymbol
    colnames(selected_counts) <- colmetadata$gene
    binary <- t(table(filtered[,c('sgene','egene')]))
    binmat <- matrix(as.numeric(as.vector(binary)), nrow=nrow(binary), ncol=ncol(binary), byrow=TRUE)
    binmat[binmat > 1] <- 1
    rownames(binmat) <- rownames(binary)
    colnames(binmat) <- colnames(binary)
    
    p <- pheatmap::pheatmap(binmat)
    print(p)

    # mnem
    k <- 3
    mnem_res <- mnem::mnem(binmat, k=k, starts=10)
    p <- plot(mnem_res)
    print(p)
	return(diffexp)
}

#' List projects we have data for
#'
#' @return List containing data objects defined by this package
#' @export
listprojects <- function() {
    projects <- c("datlinger2017_stim", "datlinger2017_unstim", "dixit2016_K562_highmoi", "dixit2016_K562_lowmoi", "jackson2019_CSTARVE", "jackson2019_DIAUXY", "jackson2019_MMD", "jackson2019_MMEtOH", "jackson2019_NLIM-GLN", "jackson2019_NLIM-NH4", "jackson2019_NLIM-PRO", "jackson2019_NLIM-UREA", "jackson2019_RAPA", "jackson2019_YPD", "jackson2019_YPEtOH", "jaitin2016")
    return(projects)
}

#' Do quality control on a SingleCellExperiment
#'
#' @param ps A SingleCellExperiment object
#' @return SingleCellExperiment object after quality control trimming
#' @export
qualitycontrol <- function(ps, plots=TRUE) {
    # get rid of rows that have never been observed
    keep_feature <- rowSums(DESeq2::counts(ps)) > 0
    summary(keep_feature) # all features are observed
    # but some are NA's
    ps <- ps[! is.na(keep_feature), ]

    # tag mitochondrial genes
    mt_genes <- grepl("MT", SummarizedExperiment::rowData(ps)$chromosome_name)
    ps <- scater::calculateQCMetrics(ps, feature_controls = list(MT = mt_genes), percent_top = 100)
    colnames(SummarizedExperiment::colData(ps))

    # experiment had no spike ins, so don't need to mark those

    if(plots) {
        p1 <- ggplot2::ggplot(as.data.frame(SummarizedExperiment::colData(ps)), ggplot2::aes(x = log10_total_counts_endogenous, y = 1, fill = 1)) +
            ggridges::geom_density_ridges(scale = 4) + 
            ggplot2::ggtitle(expression(log[10]*" library size")) + 
            ggplot2::xlab("") + ggplot2::ylab("") +
            ggplot2::scale_y_discrete(expand = c(0.01, 0)) + # will generally have to set the `expand` option
            ggridges::theme_ridges() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, face=1)) 

        # total detected genes
        p2 <- ggplot2::ggplot(as.data.frame(SummarizedExperiment::colData(ps)), ggplot2::aes(x = total_features_by_counts, y = 1,
                                                fill = 1)) +
            ggridges::geom_density_ridges(scale = 4) + 
            ggplot2::ggtitle("total features by counts detected") +
            ggplot2::xlab("") + ggplot2::ylab("") +
            ggplot2::scale_y_discrete(expand = c(0.01, 0)) +
            ggridges::theme_ridges() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, face=1))

        # percentage reads in mitochondrial genes
        p3 <- ggplot2::ggplot(as.data.frame(SummarizedExperiment::colData(ps)), ggplot2::aes(x = pct_counts_MT, y = 1,
                                                fill = 1)) +
            ggridges::geom_density_ridges(scale = 4) +
            ggplot2::ggtitle("% reads in mitochondrial genes") +
            ggplot2::xlab("") + ggplot2::ylab("") +
            ggplot2::scale_y_discrete(expand = c(0.01, 0)) +
            ggridges::theme_ridges() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, face=1))

        # percentage of counts taken by the top 100 expressed features
        p4 <- ggplot2::ggplot(as.data.frame(SummarizedExperiment::colData(ps)), ggplot2::aes(x = pct_counts_in_top_100_features, y = 1,
                                                fill = 1)) +
            ggridges::geom_density_ridges(scale = 4) + 
            ggplot2::ggtitle("% counts from top 100 expressed genes") +
            ggplot2::xlab("") + ggplot2::ylab("") +
            ggplot2::scale_y_discrete(expand = c(0.01, 0)) +
            ggridges::theme_ridges() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, face=1))

        print(scater::multiplot(plotlist = list(p1, p2, p3, p4), cols = 2))
    }

    # QC data by removing bad quality cells by removing outliers
    # the data looks like it might have already been qc'ed?
    bad_qual_mt <- scater::isOutlier(ps$pct_counts_MT, nmads=3, type="both")
    bad_qual_tf <- scater::isOutlier(ps$total_features_by_counts, nmads=3, type="both")
    bad_qual_tc <- scater::isOutlier(ps$log10_total_counts_endogenous, nmads=3, type="both")
    bad_qual_pc <- scater::isOutlier(ps$pct_counts_in_top_100_features, nmads=3, type="both")
    ps_GQ <- ps[,!(bad_qual_mt | bad_qual_tf | bad_qual_tc | bad_qual_pc)]
    # do we lose any features because cells have been removed?
    ps_GQ <- ps_GQ[rowMeans(DESeq2::counts(ps_GQ))>0,] # no op

    return(ps_GQ)

}

#' Calculate differential expression using Seurat
#'
#' @param countsdata Sparse matrix containing expression data
#' @param rowmetadata Dataframe describing the rows of countsdata (reporter genes)
#' @param colmetadata Dataframe describing the columns of countsdata (perturbed genes)
#' @return Dataframe containing the differential expression
#' @export
diffexp_seurat <- function(countsdata, rowmetadata, colmetadata) {
    perturbed.genes <- unique(colmetadata$gene)
    diffexplist <- list()
    for (gene in as.character(perturbed.genes)) {
        if (startsWith(gene, "CTRL")) {
            next
        }
        creb1 <- which(startsWith(as.character(colmetadata$gene), gene))
        ctrl <- which(startsWith(as.character(colmetadata$gene), "CTRL"))

        ctrl_counts <- countsdata[,ctrl]
        ctrl_colmetadata <- colmetadata[ctrl,]
        creb1_counts <- countsdata[,creb1]
        creb1_colmetadata <- colmetadata[creb1,]

        colnames(ctrl_counts) <- ctrl_colmetadata$cell
        colnames(creb1_counts) <- creb1_colmetadata$cell
        rownames(ctrl_counts) <- rowmetadata$origsymbol
        rownames(creb1_counts) <- rowmetadata$origsymbol

        # Set up control object
        ctrl <- Seurat::CreateSeuratObject(counts = ctrl_counts, project = "CTRL", min.cells=5)
        ctrl$stim <- "CTRL"
        ctrl <- subset(x = ctrl, subset = nFeature_RNA > 50)
        ctrl <- Seurat::NormalizeData(object = ctrl, verbose = FALSE)
        ctrl <- Seurat::FindVariableFeatures(object = ctrl, selection.method = "vst", nfeatures = 5000)

        # Set up stimulated object
        stim <- Seurat::CreateSeuratObject(counts = creb1_counts, project = "STIM", min.cells=5)
        stim$stim <- "STIM"
        stim <- subset(x = stim, subset = nFeature_RNA > 50)
        stim <- Seurat::NormalizeData(object = stim, verbose = FALSE)
        stim <- Seurat::FindVariableFeatures(object = stim, selection.method = "vst", nfeatures = 5000)

        # method 1

        #merged <- merge(ctrl, y=stim, add.cell.ids = c("CTRL", "STIM"), project="PS")
        #diffexpgenesmerged <- Seurat::FindMarkers(merged, ident.1 = "STIM", ident.2 = "CTRL" , lfc.threshold=0.1, test.use="MAST")

        # method 2

        anchors <- Seurat::FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
        combined <- Seurat::IntegrateData(anchorset = anchors, dims = 1:20)
        Seurat::DefaultAssay(object = combined) <- "integrated"

        # Run the standard workflow for visualization and clustering
        combined <- Seurat::ScaleData(object = combined, verbose = FALSE)
        combined <- Seurat::RunPCA(object = combined, npcs = 20, verbose = FALSE)
        # t-SNE and Clustering
        combined <- Seurat::RunTSNE(object = combined, reduction = "pca", dims = 1:20)
        combined <- Seurat::FindNeighbors(object = combined, reduction = "pca", dims = 1:20)
        combined <- Seurat::FindClusters(combined, resolution = 0.5)
        Seurat::DimPlot(combined, reduction="tsne", split.by="stim")

        combined$celltype.stim <- paste(Seurat::Idents(combined), combined$stim, sep = "_")
        combined$celltype <- Seurat::Idents(combined)
        nclusters <- length(levels(combined$celltype))
        Seurat::Idents(combined) <- "celltype.stim"

        allmarkers <- Seurat::FindAllMarkers(combined)
        allmarkers$sgene <- gene
        allmarkers$egene <- allmarkers$gene
        rownames(allmarkers) <- NULL
        diffexplist[[gene]] <- allmarkers

    }
    
    de <- do.call(rbind, diffexplist)
    rownames(de) <- NULL

    return(de)
}

#' Generate heatmap of knockdown efficiencies
#'
#' @param countsdata Sparse matrix containing expression data
#' @param rowmetadata Dataframe describing the rows of countsdata (reporter genes)
#' @param colmetadata Dataframe describing the columns of countsdata (perturbed genes)
#' @return pheatmap object
#' @export
knockdown_efficiency <- function() {
    data(colmetadata.dixit2016_K562_lowmoi)
    data(rowmetadata.dixit2016_K562_lowmoi)
    data(counts.dixit2016_K562_lowmoi)
    colmetadata <- colmetadata.dixit2016_K562_lowmoi
    rowmetadata <- rowmetadata.dixit2016_K562_lowmoi
    countsdata <- counts.dixit2016_K562_lowmoi

    mat <- as.matrix(countsdata)
    controls <- c("INTERGENIC1216445", "INTERGENIC1144056", "INTERGENIC216151")
    controlcells <- which(colmetadata$perturbation %in% controls)
    guides <- sort(unique(colmetadata$perturbation))
    sgenes <- sort(unique(colmetadata[!colmetadata$perturbation %in% controls,]$gene))

    pvalues <- matrix (nrow=length(sgenes), ncol=length(guides))
    rownames(pvalues) <- sgenes
    colnames(pvalues) <- guides
    
    percell <- matrix(nrow=length(sgenes), ncol=ncol(countsdata))
    rownames(percell) <- sgenes
    colnames(percell) <- rownames(colmetadata)

    for (egene in sgenes) {
        for (guide in guides) {
            #cells in the controls form the distribution
            if (egene %in% rownames(mat)) {
                controldist <- mat[egene, controlcells]
                guidecells <- which(colmetadata$perturbation %in% guide)
                guidedist <- mat[egene, guidecells]
                pval <- t.test(controldist, guidedist)$p.value
                pvalues[egene, guide] <- pval
                print(paste(egene, guide, pval, sep=" "))
            } else {
                pvalues[egene, guide] <- 1
                print(paste(egene, "not found in mat", sep=" "))
            }
        }
    }
    fraction <- 0.2
    custompalette <- rev(c(colorRampPalette(c("#ece7f2", "#2b8cbe"))(100*fraction), rep("#2b8cbe", 100*(1-fraction))))
    p <- pheatmap::pheatmap(log(pvalues, 10), cluster_cols=FALSE, cluster_rows=FALSE, color = custompalette)
    print(p)
    return(p)
}

#' Perform differential expression with deseq2
#'
#' @param countsdata Sparse matrix containing expression data
#' @param rowmetadata Dataframe describing the rows of countsdata (reporter genes)
#' @param colmetadata Dataframe describing the columns of countsdata (perturbed genes)
#' @return list of differentially expressed genes
#' @export
diffexp_deseq <- function(countsdata, rowmetadata, colmetadata) {
    controls <- c("INTERGENIC1216445", "INTERGENIC1144056", "INTERGENIC216151")
    selected.genes <- unique(colmetadata[!colmetadata$perturbation %in% controls,]$perturbation)
    selected.genes <- factor(selected.genes)

    coldata <- data.frame(row.names=rownames(colmetadata), colmetadata[,c("perturbation", "gene")])
    
    keepcols <- sample(1:nrow(coldata), 100)
    coldata <- coldata[keepcols,]

    colnames(coldata) <- c("perturbation", "gene")
    coldata$perturbation <- as.character(coldata$perturbation)
    coldata[coldata$perturbation %in% controls,]$perturbation <- "CONTROL"
    coldata$gene <- NULL
    coldata$perturbation <- factor(coldata$perturbation)

    perturbation <- coldata$perturbation
    countsdata <- floor(countsdata*1000)
    countsdata <- countsdata[,keepcols]

    dds <- DESeq2::DESeqDataSetFromMatrix(countData=countsdata, 
                                  colData=coldata, 
                                  design=~perturbation
                                  )
    
    # just pulling out the normalized counts (+divided by sizeFactor)
    dds <- DESeq2::estimateSizeFactors(dds)
    DESeq.norm.counts <- DESeq2::counts(dds, normalized=TRUE)
    
    # DGE analysis
    dds <- DESeq2::DESeq(dds)

    # to make pairwise comparisons
    # has to be done by looping 
    # over the contrasts (pairwise)
    res.list <- list()
    for(i in 1:length(selected.genes)){
        if (i %in% DESeq2::resultsNames(dds)) {
            res.list[[i]] <- DESeq2::results(dds, 
                alpha = 0.05, # alpha refers to FDR cutoff
                contrast = c("perturbation",paste0(selected.genes[i]),"CONTROL")
            ) 
        }
    } 
    names(res.list) <- selected.genes

    shrink.list <- list()
    for(i in 1:length(selected.genes)){  
        shrink.list[[i]] <- DESeq2::lfcShrink(dds, 
            contrast = c("perturbation",paste0(selected.genes[i]),"CONTROL"),
            res = res.list[[i]] # this is important, because otherwise the settings above will be neglected - e.g. the p.adjust <0.05 setting
        )
    }
    names(shrink.list) <- selected.genes

    guides <- selected.genes
    genes <- unlist(lapply(gsub("sg", "", selected.genes), function(x) {strsplit(x, "_")[[1]][[1]]}))

    for (guide in selected.genes) {
        gene <- strsplit(gsub("sg", "", guide), "_")[[1]][[1]]
        lfc <- shrink.list[[guide]][gene,]$log2FoldChange
        padj <- shrink.list[[guide]][gene,]$padj
        print(guide)
        print(padj)
    }

    return(shrink.list)
}
