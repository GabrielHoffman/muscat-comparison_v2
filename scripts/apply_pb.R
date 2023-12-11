suppressMessages({
    library(limma)
    library(muscat)
    library(scater)
    library(sctransform)
    library(SingleCellExperiment)
    library(dreamlet)
    library(edgeR)
    library(matrixStats)
    library(dreamlet)
    library(tidyverse)
    library(Matrix)
})

apply_pb <- function(sce, pars, ds_only = TRUE) {
    t1 <- system.time({
        a <- pars$assay
        if (!ds_only) 
            assay(sce, a) <- switch(a, 
                counts = counts(sce),
                cpm = calculateCPM(counts(sce)),
                logcounts = normalizeCounts(computeLibraryFactors(sce)),
                vstresiduals = vst(counts(sce), show_progress = FALSE)$y)
        pb <- aggregateData(sce, a, fun = pars$fun, scale = pars$scale)
    })[[3]]
    t2 <- system.time({
        if( pars$method %in% c("dreamlet_delta", "dreamlet_ncells", "dreamlet_none") ){

            # use dreamlet pseudobulk command here
            pb <- aggregateToPseudoBulk(sce, a, fun = pars$fun, scale = pars$scale, cluster_id = "cluster_id", sample_id = "sample_id")

            # Precision weights
            pc = 2
            W.list <- switch(pars$method, 
                    "dreamlet_delta" = pbWeights( sce, 
                                    sample_id = "sample_id", 
                                    cluster_id = "cluster_id", 
                                    method = "delta", 
                                    prior.count = pc), 
                    "dreamlet_ncells" = pbWeights( sce, 
                                    sample_id = "sample_id", 
                                    cluster_id = "cluster_id", 
                                    method = "ncells"), 
                    "dreamlet_none" = {w = pbWeights( sce, 
                                    sample_id = "sample_id", 
                                    cluster_id = "cluster_id", 
                                    method = "ncells");
                        lapply(w, function(x){x[] = 1; x})})

            vobj <- processAssays(pb, ~ group_id, verbose=FALSE, weightsList = W.list, min.cells=10, prior.count = pc)
            fit <- dreamlet(vobj, ~ group_id, verbose=FALSE )
            tab <- topTable(fit, coef='group_idB', number=Inf, sort.by="none")

            tab2 = with(tab, data.frame(gene = ID, cluster_id = assay, logFC, AveExpr, t, p_val=P.Value, B, contrast='B'))

            tab2$p_adj.glb = p.adjust(tab2$p_val, "BH")
            tab2$p_adj.loc = rep(NA, nrow(tab2))

            for( CT in unique(tab2$cluster_id) ){
                idx = which(tab2$cluster_id==CT)
                tab2$p_adj.loc[idx] = p.adjust(tab2$p_val[idx], "BH")
            }

            res = tab2    
        }else{

            # In order to keep the same genes for muscat as dreamlet
            # get gene/cluster pairs that are retained
            pb.tmp <- dreamlet::aggregateToPseudoBulk(sce, "counts", cluster_id = "cluster_id",sample_id = "sample_id")
            vobj <- dreamlet::processAssays(pb.tmp, ~ group_id, verbose=FALSE, min.cells=10, prior.count = 2)
            fit <- dreamlet(vobj, ~ group_id, verbose=FALSE )
            tab <- topTable(fit, coef='group_idB', number=Inf, sort.by="none")
            tab$key = with(tab, paste(assay, ID))

            res <- tryCatch(
                do.call(pbDS, c(
                    list(pb = pb, filter = "none", verbose = TRUE, min_cells=10),
                    pars[names(pars) %in% names(formals(pbDS))])),
                error = function(e) e)

            if (!inherits(res, "error")){
                 res <- dplyr::bind_rows(res$table[[1]])
            }

            # retain only gene/cluster pairs from dreamlet
            keep = with(res, paste(cluster_id, gene)) %in% tab$key
            res = res[keep,]
        }
    })[[3]]
    list(rt = c(t1, t2), tbl = res)
}




# Compared weighted variance along *rows* of X
# where each element of X has its own positive weight
rowWeightedVarsMatrix = function(X, W){

    stopifnot(dim(X) == dim(W))

    if( ! is.matrix(X) ) X <- as.matrix(X)
    if( ! is.matrix(W) ) W <- as.matrix(W)

    # scale weights to have a mean of 1
    Ws <- W / rowMeans2(W, useNames=FALSE)

    # weighted mean
    mu <- rowMeans2(X * Ws, useNames=FALSE)

    # weighted deviation
    A <- sqrt(Ws)*(X - mu)

    # sum of squares, divided by n-1
    rowSums2(A ** 2, useNames=FALSE) / (ncol(X)-1)
}


#' @importFrom dplyr tibble bind_rows
.getVarFromCounts <- function(countMatrix, lib.size, prior.count = .25){

    stopifnot( ncol(countMatrix) == length(lib.size))

    # if lib.size is NA, set it to zero
    lib.size[is.na(lib.size)] = 0

    countMatrix <- countMatrix + prior.count

    # pseudobulk
    count.gene <- rowSums2(countMatrix, useNames=FALSE)

    # normalize counts by library size
    # add pseudocount to counts here
    normCounts <- scale(countMatrix, 
                    scale = lib.size + 1, 
                    center = FALSE)

    # compute variance for each row
    # sigSq.hat <- rowVars(normCounts, useNames=FALSE)

    # weighted variance
    sigSq.hat <- rowWeightedVarsMatrix(normCounts, countMatrix)
    sigSq.hat[is.na(sigSq.hat)] <- 0

    # return values to compute variance later
    tibble(Gene = rownames(countMatrix), 
            count.gene = count.gene,
            sigSq.hat = sigSq.hat, 
            zeta = mean(lib.size^2), 
            ncell = ncol(countMatrix))
}

assignInNamespace("getVarFromCounts", .getVarFromCounts, ns="dreamlet")




